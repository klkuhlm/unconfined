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
       ! Theis (fully penetrating, confined)
       fp(1:np,1:nz) = theis(a,lap%p,nz)

    case(1)
       ! Hantush (partially penetrating, confined)
       fp(1:np,1:nz) = hantush(a,s%zD,s,lap%p,f,w)

    case(2)
       ! Boulton/Herrera (uconfined, non-physical boundary)
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

       where(spread(abs(eta) < MAXEXP ,2,nz))
          ! naive implementation of formula
          fp(1:np,1:nz) = udf*(1.0_EP - cosh(eta .X. s%zD)/ &
               & spread((1.0_EP + f%beta*eta*xi)*cosh(eta) + xi*sinh(eta),2,nz))
       elsewhere
          ! substitute cosh()&sinh() -> exp() and re-arrange 
          ! only valid when cosh(eta) = 0.5*exp(eta) numerically 
          fp(1:np,1:nz) = udf*(1.0_EP - exp(eta .X. (s%zD - 1.0))/ &
               & spread(1.0_EP + f%betaD*eta*xi + xi,2,nz))
       end where
       deallocate(eta,xi)

    case(5)
       ! Malama 2011 partial penetrating model (Neuman 74 when beta=0)
       allocate(eta(np),xi(np),udp(np,nz+1))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
       udp(1:np,1:nz+1) = hantush(a,[s%zD,1.0],s,lap%p,f,w)

       where(spread(abs(eta)<MAXEXP,2,nz))
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & cosh(eta .X. s%zD)/spread((1.0_EP + f%beta*eta*xi)*cosh(eta) + xi*sinh(eta),2,nz)
       elsewhere
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & exp(eta .X. (s%zD-1.0))/spread(1.0_EP + f%betaD*eta*xi + xi,2,nz)
       end where
       deallocate(eta,xi,udp)

    case(6)
       ! Mishra/Neuman  2010 model
       fp(1:np,1:nz) = mishraNeuman2010(a,s%zD,s,lap%p,f,w)
       
    end select

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
    integer :: np,nz
    integer, allocatable :: zLay(:)

    nz = size(zD)
    np = size(p)

    allocate(zLay(nz))
    if(size(s%zLay) == nz) then
       zLay = s%zLay
    else
       zLay = [s%zLay, 3] ! last point is at the water table (always layer 3)
    end if

    eta(1:np) = sqrt((p(:) + a**2)/f%kappa)
    xi(1:np) = eta(:)*f%alphaD/p(:)

    where(spread(zLay == 3,1,np))
       ! above well screen
       g(1,1:np,1:nz) = cosh(eta(1:np) .X. (1.0 - w%dD - zD(1:nz)))
    end where
    
    where(spread(zLay == 1,1,np))
       ! below well screen
       g(3,1:np,1:nz) = spread(exp(-eta*(1.0 - w%lD)) - &
            & (sinh(eta*w%dD) + exp(-eta)*sinh(eta*(1.0 - w%lD)))/sinh(eta),2,nz)
    end where

    g(2,1:np,1:nz) = (spread(sinh(eta*w%dD),2,nz)*cosh(eta .X. zD) + &
                    & spread(sinh(eta*(1.0-w%lD)),2,nz)*cosh(eta .X. (1.0-zD)))/&
                    & spread(sinh(eta),2,nz)

    where(spread(zLay(1:nz) == 1,1,np))
       ! below well screen (0 <= zD <= 1-lD)
       udp(1:np,1:nz) =  g(3,:,:)*cosh(eta .X. zD)
    elsewhere
       where(spread(zLay(1:nz) == 2,1,np))
          ! next to well screen (1-lD < zD < 1-dD)
          udp(1:np,1:nz) = (1.0_EP - g(2,:,:))
       elsewhere
          ! above well screen (1-lD <= zD <= 1)
          udp(1:np,1:nz) = (g(1,:,:) - g(2,:,:))
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

    ! only used in asymptotic expansion
    real(EP), dimension(size(p)) :: nuep
    complex(EP), dimension(size(p)) :: phiep

    real(DP) :: nu, B2
    real(EP), dimension(size(p)) :: delta2
    real(DP), dimension(0:3) :: beta

    complex(DP), dimension(size(p)) :: phi
    complex(EP) :: arg1
    complex(EP), dimension(size(p)) :: arg2, B1, delta1

    complex(EP), dimension(2,2,size(p)) :: aa
    complex(EP), dimension(size(p)) :: eta
    complex(EP), dimension(size(p),2) :: J,Y
    complex(DP), dimension(2) :: tmp
    
    
    np = size(p)
    nz = size(zD)

    beta(0) = f%acD*f%sigma
    beta(1) = f%lambdaD
    beta(2) = f%ak*f%b1
    beta(3) = f%akD

    eta(1:np) = sqrt((a**2 + p(1:np))/f%kappa)
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0],s,p,f,w)
    
    B1(1:np) = p(1:np)*beta(0)/f%kappa*exp(-beta(2))
    B2 = real(a**2/f%kappa,DP)

    phi(1:np) = cmplx(2.0*EYE*sqrt(B1/beta(1)**2)*exp(beta(1)*f%usLD/2.0),kind=DP)
    nu = sqrt((beta(3)**2 + 4.0*B2)/beta(2)**2)

    phiep(1:np) = 2.0_EP*EYE*sqrt(B1/beta(1)**2)*exp(beta(1)*f%usLD/2.0_EP)
    nuep(1:2) = real([nu,nu+1.0],EP)

    call cbesj(z=phi(1),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
    if (ierr == 4 .or. ierr == 2 .or. abs(tmp(2)) >= huge(1.0_DP)) then
       J(:,:) = besselJAsymptoticOrder(nuep,phiep)
    else
       do i= 1,np
          call cbesj(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
          if (ierr > 0 .and. ierr /= 3) then
             print *, 'ERROR: CBESJ (zD=LD)',phi(i),nu,i,ierr,nzero
             stop
          else
             J(i,1:2) = tmp(1:2)
          end if
       end do
    end if

    call cbesy(z=phi(1),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
    if (ierr == 4 .or. ierr == 2 .or. abs(tmp(2)) >= huge(1.0_DP)) then
       Y(:,:) = besselYAsymptoticOrder(nuep,phiep)
    else
       do i=1,np
          call cbesy(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
          if (ierr > 0 .and. ierr /= 3) then
             print *, 'ERROR: CBESY (zD=LD)',phi(i),nu,i,ierr,nzero
             stop
          else
             Y(i,1:2) = tmp(1:2)
          end if
       end do
    end if
    
    arg1 = real(beta(3),EP) + real(nu*beta(1),EP)
    arg2(1:np) = beta(1)*phi(1:np)
    aa(1,1,1:np) = arg1*J(:,1) - arg2(:)*J(:,2)
    aa(1,2,1:np) = arg1*Y(:,1) - arg2(:)*Y(:,2)
    phi(1:np) = cmplx(2.0*EYE*sqrt(B1/beta(1)**2),kind=DP)


    call cbesj(z=phi(1),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
    if (ierr == 4 .or. ierr == 2 .or. abs(tmp(2)) >= huge(1.0_DP)) then
       J(:,:) = besselJAsymptoticOrder(nuep,phiep)
    else
       do i= 1,np
          call cbesj(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
          if (ierr > 0 .and. ierr /= 3) then
             print *, 'ERROR: CBESJ (zD=LD)',phi(i),nu,i,ierr,nzero
             stop
          else
             J(i,1:2) = tmp(1:2)
          end if
       end do
    end if

    call cbesy(z=phi(1),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
    if (ierr == 4 .or. ierr == 2 .or. abs(tmp(2)) >= huge(1.0_DP)) then
       Y(:,:) = besselYAsymptoticOrder(nuep,phiep)
    else
       do i=1,np
          call cbesy(z=phi(i),fnu=nu,kode=2,n=2,cy=tmp(1:2),nz=nzero,ierr=ierr)
          if (ierr > 0 .and. ierr /= 3) then
             print *, 'ERROR: CBESY (zD=LD)',phi(i),nu,i,ierr,nzero
             stop
          else
             Y(i,1:2) = tmp(1:2)
          end if
       end do
    end if
    
    arg2(1:np) = beta(1)*phi(1:np)
    aa(2,1,1:np) = arg1*J(:,1) - arg2(:)*J(:,2)
    aa(2,2,1:np) = arg1*Y(:,1) - arg2(:)*Y(:,2)

    delta2(1:np) = abs(aa(1,1,:)*aa(2,2,:) - aa(1,2,:)*aa(2,1,:))
    delta1(1:np) = (aa(1,1,:)*Y(:,1) - aa(1,2,:)*J(:,1))*2.0*eta(:)* &
         & sinh(eta(:))/delta2(:) - cosh(eta(:))

    sU(1:np,1:nz) = spread(sH(1:np,nz+1)/delta1(:),2,nz)*cosh(eta(:) .X. s%zD(:))
    sD(1:np,1:nz) = sH(:,1:nz) + sU(:,:)

    write(*,999) ' nu:',nu,' z:',phi(1:3),' sU:',su(1:3,1),' D1:',delta1(1:3),&
         & ' D2:',delta2(1:3),' sH:',sH(1:3,1)

999 format(A,ES11.3E3,3(A,3('(',ES12.3E4,',',ES12.3E4,')')),A,3(ES12.3E4,1X),A,&
         & 3('(',ES12.3E4,',',ES12.3E4,')'))

  end function mishraNeuman2010

  function besselJAsymptoticOrder(nu,z) result(J)
    use constants, only : EP, PIEP, E, SQRT2
    real(EP), dimension(2), intent(in) :: nu ! order
    complex(EP), dimension(:), intent(in) :: z ! argument
    complex(EP), dimension(size(z),2) :: J
    integer :: i,k

    ! evaluate the asymptotic expansion using extended range
    ! so it can be used when the double-precision Amos routines
    ! return overflow or underflow

    ! http://dlmf.nist.gov/10.19#i
    forall (k=1:size(z), i=1:2)
       J(k,i) = 1.0/sqrt(2.0*PIEP*nu(i))*(E*z(k)/(2.0*nu(i)))**nu(i)
    end forall

  end function  besselJAsymptoticOrder

  function besselYAsymptoticOrder(nu,z) result (Y)
    use constants, only : EP, PIEP, E, SQRT2
    real(EP), dimension(2), intent(in) :: nu ! order
    complex(EP), dimension(:), intent(in) :: z ! argument
    complex(EP), dimension(size(z),2) :: Y
    integer :: i,k

    ! evaluate the asymptotic expansion using extended range
    ! so it can be used when the double-precision Amos routines
    ! return overflow or underflow

    ! http://dlmf.nist.gov/10.19#i
    forall (k=1:size(z), i=1:2)
       Y(k,i) =      -SQRT2/(PIEP*nu(i))*(E*z(k)/(2.0*nu(i)))**(-nu(i))
    end forall

  end function besselYAsymptoticOrder


end module laplace_hankel_solutions




