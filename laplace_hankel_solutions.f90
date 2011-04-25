 module laplace_hankel_solutions
  implicit none

  private
  public :: lap_hank_soln
contains

  function lap_hank_soln(a,rD,np,nz,w,f,s,lap) result(fp)
    use constants, only : DP, EP, MAXEXP
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.), cosh, sinh 

    implicit none
    
    real(EP), intent(in) :: a
    real(DP), intent(in) :: rD
    integer, intent(in) :: np,nz
    type(invLaplace), intent(in) :: lap
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(np,nz) :: fp

    complex(EP), allocatable :: eta(:), xi(:), udp(:,:), udf(:,:)

    intrinsic :: bessel_j0

    select case(s%model)
    case(0)
       ! Theis (fully penetrating confined)
       fp(1:np,1:nz) = theis(a,lap%p,nz)

    case(1)
       ! Hantush (partially penetrating confined)
       fp(1:np,1:nz) = hantush(a,s%zD,s,lap%p,f,w)

    case(2)
       ! Boulton/Herrera (uconfined, non-physical boundary condition)
       stop 'ERROR Boulton/Herrera model not implemented yet'

    case(3)
       ! Moench (unconfined w/ wellbore storage)
       stop 'ERROR Moench model not implmented yet'

    case(4)
       ! Malama 2011 fully penetrating model (Neuman 72 when beta=0)
       allocate(eta(np),xi(np),udf(np,nz))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
       udf(1:np,1:nz) = theis(a,lap%p,nz)

!!$       where(spread(abs(eta) < MAXEXP ,2,nz))
          ! naive implementation of formula
          fp(1:np,1:nz) = udf*(1.0_EP - cosh(eta .X. s%zD)/ &
               & spread((1.0_EP + f%beta*eta*xi)*cosh(eta) + xi*sinh(eta),2,nz))
!!$       elsewhere
!!$          ! substitute cosh()&sinh() -> exp() and re-arrange 
!!$          ! only valid when cosh(eta) = 0.5*exp(eta) numerically 
!!$          fp(1:np,1:nz) = udf*(1.0_EP - exp(eta .X. (s%zD - 1.0))/ &
!!$               & spread(1.0_EP + f%betaD*eta*xi + xi,2,nz))
!!$       end where
       deallocate(eta,xi)

    case(5)
       ! Malama 2011 partial penetrating model (Neuman 74 when beta=0)
       allocate(eta(np),xi(np),udp(np,nz+1))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
       udp(1:np,1:nz+1) = hantush(a,[s%zD,1.0],s,lap%p,f,w)

!!$       where(spread(abs(eta)<MAXEXP,2,nz))
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & cosh(eta .X. s%zD)/spread((1.0_EP + f%beta*eta*xi)*cosh(eta) + xi*sinh(eta),2,nz)
!!$       elsewhere
!!$          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
!!$               & exp(eta .X. (s%zD-1.0))/spread(1.0_EP + f%betaD*eta*xi + xi,2,nz)
!!$       end where
       deallocate(eta,xi,udp)

    case(6)
       ! Mishra/Neuman  2010 model
       fp(1:np,1:nz) = mishraNeuman2010(a,s%zD,s,lap%p,f,w)
       
    end select

    ! apply common Hankel-transform and time-behavior factors
    fp(1:np,1:nz) = a*bessel_j0(a*rD)*fp(:,:)*spread(lapTime(lap),2,nz)

  end function lap_hank_soln
  
  function theis(a,p,nz) result(udf)
    use constants, only : EP
    real(EP), intent(in) :: a
    complex(EP), dimension(:), intent(in) :: p
    integer, intent(in) :: nz
    complex(EP), dimension(size(p),nz) :: udf

    udf(:,:) = spread(2.0_EP/(p(:) + a**2),dim=2,ncopies=nz)

  end function theis

  function hantush(a,zD,s,p,f,w) result(udp)
    ! implemented in the form given in Malama,Kuhlman & Barrash 2008 
    use constants, only : DP, EP, MAXEXP
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.), cosh, sinh
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(size(p),size(zd)) :: udp
    complex(EP), dimension(size(p)) :: eta, xi
    complex(EP), dimension(3,size(p),size(zd)) :: g
    complex(EP), dimension(2,size(p),size(zd)) :: ff
    integer :: np,nz
    integer, allocatable :: zLay(:,:)

    nz = size(zD) ! # z vals requested 
    np = size(p)

    allocate(zLay(np,nz))

    ! are # zvals requested same as given in program input?
    if(size(s%zLay) == nz) then
       zLay(1:np,1:nz) = spread(s%zLay,1,np)
    else
       ! last point is at the water table for boundary conditions
       ! (top boundary condition is always layer 3)
       zLay(1:np,1:nz) = spread([s%zLay, 3],1,np)
    end if

    eta(1:np) = sqrt((p(:) + a**2)/f%kappa)
    xi(1:np) = eta(:)*f%alphaD/p(:)

    ff(1,1:np,1:nz) = spread(sinh(eta*w%dD),dim=2,ncopies=nz)
    ff(2,1:np,1:nz) = spread(sinh(eta*(1.0-w%lD)),dim=2,ncopies=nz)

    where(zLay == 3)
       ! above well screen
       g(1,1:np,1:nz) = cosh(eta(1:np) .X. (1.0 - w%dD - zD(1:nz)))
    end where

    where(zLay == 1 .or. zLay == 2)
       g(2,1:np,1:nz) = (ff(1,:,:)*cosh(eta .X. zD) + ff(2,:,:)*cosh(eta .X. (1.0-zD)))/&
            & spread(sinh(eta),2,nz)
    end where
    
    where(zLay == 1)
       ! below well screen
       g(3,1:np,1:nz) = spread(exp(-eta*(1.0 - w%lD)) - &
            & (ff(1,:,1) + exp(-eta)*ff(2,:,1))/sinh(eta),2,nz)
    end where

    where(zLay == 1)
       ! below well screen (0 <= zD <= 1-lD)
       udp(1:np,1:nz) =  g(3,:,:)*cosh(eta .X. zD)
    elsewhere
       where(zLay == 2)
          ! next to well screen (1-lD < zD < 1-dD)
          udp(1:np,1:nz) = 1.0_EP - g(2,:,:)
       elsewhere
          ! zLay == 3
          ! above well screen (1-lD <= zD <= 1)
          udp(1:np,1:nz) = g(1,:,:) - g(2,:,:)
       end where
    end where

    udp(1:np,1:nz) = udp(:,:)*theis(a,p,nz)/w%bD
  end function hantush
  
  function mishraNeuman2010(a,zD,s,p,f,w) result(sD)
    use constants, only : DP, EP, EYE, PIEP, E, SQRT2
    use types, only : well, formation, solution
    use utility, only : operator(.X.), cosh, sinh
    use cbessel, only : cbesj, cbesy ! Amos routines
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f

    complex(EP), dimension(size(p),size(zD)) :: sD, sU
    complex(EP), dimension(size(p),size(zD)+1) :: sH
    integer :: i, nzero, ierr, np, nz

    ! extended precision version used in asymptotic expansion
    real(EP), dimension(size(p)) :: nuep
    complex(EP), dimension(size(p)) :: phiep

    real(DP) :: nu
    real(DP), dimension(0:3) :: beta

    complex(EP), dimension(size(p)) :: delta2
    real(EP) :: B2

    complex(DP), dimension(size(p)) :: phi
    complex(DP), dimension(2) :: tmp

    complex(EP) :: arg1
    complex(EP), dimension(size(p)) :: arg2, B1, delta1
    complex(EP), dimension(2,2,size(p)) :: aa
    complex(EP), dimension(size(p)) :: eta
    complex(EP), dimension(size(p),2) :: J,Y
    
    np = size(p)
    nz = size(zD)

    beta(0) = f%ac*f%Sy/f%Ss
    beta(1) = f%lambdaD  ! b*(ak-ac)
    beta(2) = f%ak*f%b1  ! ak*(psi_a - psi_k);
    beta(3) = f%akD

    eta(1:np) = sqrt((a**2 + p(1:np))/f%kappa)
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0],s,p,f,w)
    
    B1(1:np) = p(:)*beta(0)*exp(-beta(2))/f%kappa
    B2 = (a**2)/f%kappa

    ! compute v1
    phiep(1:np) = EYE*sqrt(4.0*B1/(beta(1)**2))*exp(beta(1)*f%usLD/2.0_EP)
    phi(1:np) = cmplx(phiep,kind=DP) 

    nuep(1) = sqrt((beta(3)**2 + 4.0*B2)/beta(1)**2)
    nuep(2) = nuep(1) + 1.0_EP
    nu = real(nuep(1),kind=DP) 
    
    do i= 1,np
       ! kode=2 is scaled BF
       call cbesj(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
       if (ierr > 0 .and. ierr /= 3) then
          print *, 'ERROR: CBESJ (zD=LD) z=',phi(i),' nu=',nu,' i,ierr,nz:',i,ierr,nzero
!!$             stop
       else
          J(i,1:2) = tmp(1:2)
       end if
       call cbesy(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
       if (ierr > 0 .and. ierr /= 3) then
          print *, 'ERROR: CBESY (zD=LD) z=',phi(i),' nu=',nu,' i,ierr,nz:',i,ierr,nzero
!!$             stop
       else
          Y(i,1:2) = tmp(1:2)
       end if
    end do
    
    ! compute v3
    arg1 = real(beta(3),EP) + nuep(1)*beta(1)
    arg2(1:np) = beta(1)*phiep(1:np)

    aa(1,1,1:np) = arg1*J(:,1) - arg2(:)*J(:,2)
    aa(1,2,1:np) = arg1*Y(:,1) - arg2(:)*Y(:,2)

    ! compute v2
    phiep(1:np) = EYE*sqrt(4.0*B1/beta(1)**2)
    phi(1:np) = cmplx(phiep,kind=DP)

    ! kode=2 is scaled BF
    do i= 1,np
       call cbesj(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
       if (ierr > 0 .and. ierr /= 3) then
          print *, 'ERROR: CBESJ (zD=0) z=',phi(i),' nu=',nu,' i,ierr,nz:',i,ierr,nzero
!!$             stop
       else
          J(i,1:2) = tmp(1:2)
       end if
       call cbesy(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
       if (ierr > 0 .and. ierr /= 3) then
          print *, 'ERROR: CBESY (zD=0) z=',phi(i),' nu=',nu,' i,ierr,nz:',i,ierr,nzero
!!$             stop
       else
          Y(i,1:2) = tmp(1:2)
       end if
    end do
    
    arg2(1:np) = beta(1)*phiep(1:np)
    aa(2,1,1:np) = arg1*J(:,1) - arg2(:)*J(:,2)
    aa(2,2,1:np) = arg1*Y(:,1) - arg2(:)*Y(:,2)

    ! product of phi(0) and phi(lD) normalizations
    delta2(1:np) = (aa(1,1,:)*aa(2,2,:) - aa(1,2,:)*aa(2,1,:)) 

    ! products of normalizations cancel
    delta1(1:np) = (aa(1,1,:)*Y(:,1) - aa(1,2,:)*J(:,1))/delta2(:)* &
         & 2.0*eta(:)*sinh(eta(:)) - cosh(eta(:))

    sU(1:np,1:nz) = spread(sH(1:np,nz+1)/delta1(1:np),2,nz)*cosh(eta(1:np) .X. s%zD(1:nz))
    sD(1:np,1:nz) = sH(1:np,1:nz) + sU(1:np,1:nz)

    if (s%quiet > 1) then
       write(*,999) ' nu:',nu,' phi:',phiep(1:3),' sU:',su(1:3,1),' D1:',delta1(1:3),&
            & ' D2:',delta2(1:3),' sH:',sH(1:3,1)
    end if

999 format(A,ES11.3E3,5(A,3('(',ES12.3E4,',',ES12.3E4,')')))
  end function mishraNeuman2010
  
!!$  subroutine besselJYAsymptoticArg(nu,z,J,Y)
!!$    use constants, only : EP, PIEP, PIOV2EP, PIOV4EP, E, SQRT2
!!$    integer, parameter :: M = 30, hM = M/2-1
!!$    real(EP), dimension(2), intent(in) :: nu ! order
!!$    complex(EP), dimension(:), intent(in) :: z ! argument
!!$    complex(EP), dimension(size(z),2) :: J, Y
!!$    integer :: i,k,nz
!!$
!!$    real(EP), dimension(0:M,2) :: a
!!$    complex(EP), dimension(size(z),2) :: omega
!!$    real(EP) :: gam, ek
!!$    integer, dimension(0:hM) :: iv ! integer vector
!!$    real(EP), dimension(0:hM,2) :: sv  ! sign vector
!!$    complex(EP), dimension(2,size(z)) :: arg1,arg2
!!$
!!$    nz = size(z)
!!$    a = 1.0_EP; gam = 1.0_EP; ek = 1.0_EP
!!$    do i=1,M
!!$       ! build up k! and 8**k as we step
!!$       gam = gam*i
!!$       ek = ek*8.0_EP
!!$       a(i:M,1:2) = a(i:M,1:2)*spread((4.0*nu(1:2) - i**2)/(gam*ek),1,M-i+1)
!!$    end do
!!$    
!!$    omega(1:nz,1:2) = spread(z(:),2,2) - spread(nu(:)*PIOV2EP,1,nz) - PIOV4EP
!!$    forall(i=0:hM) iv(i) = i
!!$    sv(0:hM,1:2) = spread(1.0_EP**iv(0:hM),2,2)
!!$
!!$    arg1(1:2,1:nz) = sum(spread(sv(0:hM,1:2)*a(2*iv,1:2),3,nz)/&
!!$         & spread(spread(z(1:nz),1,2),1,hM+1)**&
!!$         & spread(spread(2*iv,2,2),3,nz),dim=1)
!!$
!!$    arg2(1:2,1:nz) = sum(spread(sv(0:hM,1:2)*a(2*iv+1,1:2),3,nz)/&
!!$         & spread(spread(z(1:nz),1,2),1,hM+1)**&
!!$         & spread(spread(2*iv+1,2,2),3,nz),dim=1)
!!$
!!$    ! http://dlmf.nist.gov/10.17#i
!!$    forall(k=1:nz, i=1:2)
!!$       J(k,i) = sqrt(2.0/(PIEP*z(k)))*&
!!$            & (cos(omega(k,i))*arg1(i,k) - sin(omega(k,i))*arg2(i,k))
!!$       Y(k,i) = sqrt(2.0/(PIEP*z(k)))*&
!!$            & (sin(omega(k,i))*arg1(i,k) + cos(omega(k,i))*arg2(i,k))
!!$    end forall
!!$
!!$  end subroutine besselJYAsymptoticArg

!!$  subroutine besselJYAsymptoticOrder(nu,z,J,Y)
!!$    use constants, only : EP, PIEP, E, SQRT2
!!$    real(EP), dimension(2), intent(in) :: nu ! order
!!$    complex(EP), dimension(:), intent(in) :: z ! argument
!!$    complex(EP), dimension(size(z),2), intent(out) :: J, Y
!!$    integer :: i,k
!!$
!!$    ! evaluate the asymptotic expansion using extended range
!!$    ! so it can be used when the double-precision Amos routines
!!$    ! return overflow or underflow
!!$
!!$    ! http://dlmf.nist.gov/10.19#i
!!$    forall (k=1:size(z), i=1:2)
!!$       J(k,i) = 1.0/sqrt(2.0*PIEP*nu(i))*(E*z(k)/(2.0*nu(i)))**nu(i)
!!$       Y(k,i) = -SQRT2/(PIEP*nu(i))*(E*z(k)/(2.0*nu(i)))**(-nu(i))
!!$    end forall
!!$
!!$  end subroutine   besselJYAsymptoticOrder

end module laplace_hankel_solutions




