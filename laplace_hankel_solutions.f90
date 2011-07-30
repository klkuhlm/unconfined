 module laplace_hankel_solutions
  implicit none

  private
  public :: lap_hank_soln
contains

  function lap_hank_soln(a,rD,np,nz,w,f,s,lap) result(fp)
    use constants, only : DP, EP, MAXEXP
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.) 

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
       ! Hantush-style solution with wellbore storage
       fp(1:np,1:nz) = hantushstorage(a,s%zD,s,lap%p,f,w)

    case(3)
       ! Moench (unconfined w/ wellbore storage)
       stop 'ERROR Moench model not implmented yet'

    case(4:5)
       ! Malama 2011 partial penetrating model (Neuman 74 when beta=0)
       allocate(eta(np),xi(np),udp(np,nz+1))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
       if (s%model == 4) then
          udf(1:np,1:nz) = theis(a,lap%p,nz)
          udf(1:np,nz+1) = udf(1:np,nz)
       else
          udp(1:np,1:nz+1) = hantush(a,[s%zD,1.0],s,lap%p,f,w)
       end if
       
       where (spread(real(eta) < MAXEXP,2,nz))
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & cosh(eta .X. s%zD)/spread((1.0_EP + f%beta*eta*xi)*&
               & cosh(eta) + xi*sinh(eta),2,nz)          
       elsewhere
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & exp(eta .X. (s%zD - 1.0))/&
               & spread(1.0_EP + f%beta*eta*xi + xi,2,nz)
       end where
       deallocate(eta,xi,udp)

    case(6)
       ! Mishra/Neuman  2010 model
       fp(1:np,1:nz) = mishraNeuman2010b(a,s%zD,s,lap%p,f,w)
       
    end select

    ! apply common Hankel-transform and time-behavior factors
    fp(1:np,1:nz) = a*bessel_j0(a*rD)*fp*spread(lapTime(lap),2,nz)

  end function lap_hank_soln
  
  function theis(a,p,nz) result(udf)
    use constants, only : EP
    real(EP), intent(in) :: a
    complex(EP), dimension(:), intent(in) :: p
    integer, intent(in) :: nz
    complex(EP), dimension(size(p),nz) :: udf

    udf(:,:) = spread(2.0_EP/(p(:) + a**2),dim=2,ncopies=nz)

  end function theis

  function hantush(a,zD,sol,p,f,w) result(udp)
    ! implemented in the form given in Malama, Kuhlman & Barrash 2008 
    use constants, only : DP, EP, SMALLZ
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.) 
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: sol
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(size(p),size(zD)) :: udp
    complex(EP), dimension(size(p)) :: eta
    complex(EP), dimension(3,size(p),size(zD)) :: g
    complex(EP), dimension(2,size(p),size(zD)) :: ff
    integer, dimension(size(p),size(zD)) :: zLay
    integer :: np,nz

    nz = size(zD) ! # z vals requested of this function
    np = size(p)

    ! are # zvals requested same as given in program input?
    if(sol%nz == nz) then
       ! called directly
       zLay(1:np,1:nz) = spread(sol%zLay(:),1,np)
    else
       ! called from other method, requiring solution 
       ! at the water table for boundary conditions
       ! (top boundary condition is always layer 3)
       zLay(1:np,1:nz) = spread([sol%zLay(:),3],1,np)
    end if

    eta(1:np) = sqrt((p(:) + a**2)/f%kappa)
    
    ! above well screen
     g(1,1:np,1:nz) = cosh(eta(:) .X. (1.0 - w%dD - zD(1:nz)))
    ff(1,1:np,1:nz) = spread(sinh(eta*w%dD),2,nz)
    ff(2,1:np,1:nz) = spread(sinh(eta*(1.0 - w%lD)),2,nz)

    g(2,1:np,1:nz) = (ff(1,:,:)*cosh(eta .X. zD) + ff(2,:,:)*&
         & cosh(eta .X. (1.0 - zD)))/spread(sinh(eta),2,nz)

    ! below well screen
    g(3,1:np,1:nz) = spread(exp(-eta*(1.0 - w%lD)) - &
         & (ff(1,:,1) + exp(-eta)*ff(2,:,1))/sinh(eta),2,nz)

    where(zLay == 1)
       ! below well screen (0 <= zD <= 1-lD)
       udp(1:np,1:nz) = g(3,:,:)*cosh(eta .X. zD)
    elsewhere
       where (zLay == 2)
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
  
  function hantushstorage(a,zD,sol,p,f,w) result(u)
    use constants, only : DP, EP, PI, SMALLZ
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.) 
    use cbessel, only : cbesk ! Amos routine
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: sol
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(size(p),size(zD)) :: u, uDp

    complex(DP), dimension(0:1,size(p)) :: K
    complex(EP), dimension(size(p)) :: eta, xi, A0, uDf
    complex(EP), dimension(3,size(p),size(zD)) :: ff
    complex(EP), dimension(2,size(p),size(zD)) :: g

    real(DP) :: CDw, tDb
    integer :: np,nz,i

    ! things related to Amos library
    integer(4), parameter :: kode = 1, num = 2
    integer(4) :: nzero, ierr
    integer, dimension(size(p),size(zD)) :: zLay

    np = size(p)
    nz = size(zd)

    ! are # zvals requested same as given in program input?
    if(sol%nz == nz) then
       zLay(1:np,1:nz) = spread(sol%zLay(:),1,np)
    else
       zLay(1:np,1:nz) = spread([sol%zLay(:),3],1,np)
    end if

    ! dimensionless well storage coefficient
    CDw = w%rDw**2/(2.0*(w%l - w%d)*f%Ss)
    
    ! dimensionless well delay time
    tDb = PI*sol%rDwobs**2/(sol%sF*f%Ss)
    
    xi(1:np) = w%rDw*sqrt(p(:))
    eta(1:np) = sqrt((p(:) + a**2)/f%kappa)

    do i=1,np
       call cbesk(z=cmplx(xi(i),kind=DP),fnu=0.0,kode=kode,&
            & n=num,cy=K(0:1,i),nz=nzero,ierr=ierr)
       if (ierr > 0 .and. ierr /= 3) then
          print *, 'ERROR: CBESK z=',xi(i),&
               &' i,ierr,nz:',i,ierr,nzero
       end if
    end do
    
    A0(1:np) = 2.0/(p(:)*CDw*K(0,:) + xi(:)*K(1,:))
    uDf(1:np) = A0(:)/(p*(p + a**2)*(p*tDb + 1.0_EP))
        
    where (zLay == 3)
       g(1,1:np,1:nz) = cosh(eta(:) .X. (1.0 - w%dD - zd(:)))
    end where
    
    where (zLay == 1 .or. zLay == 2)
       ff(1,1:np,1:nz) = spread(sinh(eta(:)*w%dD),2,nz)
       ff(2,1:np,1:nz) = spread(sinh(eta(:)*(1.0 - w%lD)),2,nz)
    end where
    
    where (zLay == 2 .or. zLay == 3)
       g(2,1:np,1:nz) = (ff(1,:,:)*cosh(eta .X. zd) + &
            & ff(2,:,:)*cosh(eta .X. (1.0 - zd)))/&
            & spread(sinh(eta(:)),2,nz)
    end where

    where (zLay == 1)
       ff(3,1:np,1:nz) = spread(exp(-eta(:)*(1.0 - w%lD)) - &
            & (ff(1,:,1) + exp(-eta(:))*ff(2,:,1))/sinh(eta(:)),2,nz)
    end where

    where (zLay == 3) !!(spread(zD > 1-w%dD,2,nz))
       uDp(1:np,1:nz) = g(1,:,:) - g(2,:,:)
    elsewhere 
       where (zLay == 1) !!(spread(zD < (1.0-w%lD),2,nz))
          uDp(1:np,1:nz) = ff(3,:,:)*cosh(eta .X. zd)
       elsewhere
          ! zLay == 2
          uDp(1:np,1:nz) = 1.0_EP - g(2,:,:)
       end where
    end where
    
    u(1:np,1:nz) = spread(uDf(:)/w%bD,2,nz)*uDp(:,:)

  end function hantushstorage

  function mishraNeuman2010(a,zD,s,p,f,w) result(sD)
    use constants, only : DP, EP, EYE, PIEP, E, SQRT2
    use types, only : well, formation, solution
    use utility, only : operator(.X.) 
    use complex_bessel
    use complex_bessel2
    use cbessel, only : cbesj,cbesy ! Amos routine
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f

    complex(EP), dimension(size(p),size(zD)) :: sD, sU
    complex(EP), dimension(size(p),size(zD)+1) :: sH
    integer ::  np, nz

    real(EP) :: nuep, B2, nv
    complex(EP), dimension(size(p)) :: phiep

    real(DP), dimension(0:3) :: beta
    complex(EP) :: arg1
    complex(EP), dimension(size(p)) :: arg2, B1, delta1, delta2
    complex(EP), dimension(2,2,size(p)) :: aa
    complex(EP), dimension(size(p)) :: eta
    complex(EP), dimension(size(p),2) :: J,Y
    
    integer, parameter :: NPRINT = 2

    ! size integer expected by BF library
    integer(4), parameter :: kode = 2, num = 2
    integer(4) :: nzero, ierr
    complex(DP), dimension(size(p)) :: phi
    complex(EP), allocatable :: tmp2(:,:)
    complex(DP), dimension(2) :: tmp
    real(DP) :: nu
    integer :: i

    intrinsic :: isnan

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
    phiep(1:np) = EYE*sqrt(4.0*B1/(beta(1)**2))*exp(0.5_DP*beta(1)*f%usLD)
    nuep = sqrt((beta(3)**2 + 4.0*B2)/beta(1)**2)
 
    phi(1:np) = cmplx(phiep(1:np),kind=DP)
    nu = real(nuep,kind=DP)

    if (all(abs(phi(:)) < 2.0*(nu+1.0))) then
       call cjy(nuep,    phiep(:),J(:,1),Y(:,1))
       call cjy(nuep+1.0,phiep(:),J(:,2),Y(:,2))
    else
       do i= 1,np
          if (2.0*(nu+1.0) > abs(phi(i))) then
             call cjy(nuep,    phiep(i:i),J(i,1),Y(i,1))
             call cjy(nuep+1.0,phiep(i:i),J(i,2),Y(i,2))
          else
             call cbesj(z=phi(i),fnu=nu,kode=kode,n=num,cy=tmp(1:2),&
                  & nz=nzero,ierr=ierr)
             if (ierr > 0  .and. ierr /= 3 .and. s%quiet > 1) then
                print *, 'ERROR: CBESJ (zD=LD) z=',phi(i),' nu=',nu,&
                     &' i,ierr,nz:',i,ierr,nzero  

             else
                J(i,1:2) = tmp(1:2)
             end if
             call cbesy(z=phi(i),fnu=nu,kode=kode,n=num,cy=tmp(1:2),&
                  &nz=nzero,ierr=ierr)
             if (ierr > 0  .and. ierr /= 3 .and. s%quiet > 1) then
                print *, 'ERROR: CBESY (zD=LD) z=',phi(i),' nu=',nu,&
                     &' i,ierr,nz:',i,ierr,nzero
             else
                Y(i,1:2) = tmp(1:2)
             end if

             if (any(isnan(abs(J(i,1:2)))) .or. any(isnan(abs(Y(i,1:2))))) then
                allocate(tmp2(4,0:int(nu)+2))                   
                call cjyva(nuep+1.0,phiep(i),nv,tmp2(1,:),tmp2(2,:),tmp2(3,:),tmp2(4,:))
                J(i,1:2) = tmp2(1,int(nv)-1:int(nv))
                Y(i,1:2) = tmp2(3,int(nv)-1:int(nv))
                deallocate(tmp2)
                print *, 'zD=LD',nuep+1.0,phiep(i),J(i,1),Y(i,1)
             end if
          end if
       end do
    end if
    
    ! compute v3
    arg1 = real(beta(3),EP) + nuep*beta(1)
    arg2(1:np) = beta(1)*phiep(1:np)

    aa(1,1,1:np) = arg1*J(:,1) - arg2(:)*J(:,2)
    aa(1,2,1:np) = arg1*Y(:,1) - arg2(:)*Y(:,2)

    if (s%quiet > 1) then
       write(*,998) 'zD=LD phi:',phiep(1:NPRINT),'J:',J(1:NPRINT,1),&
            & 'Y:',Y(1:NPRINT,1),'a(1,1)',aa(1,1,1:NPRINT),'a(1,2)',aa(1,2,1:NPRINT)
    end if

998 format(5(A,2('(',ES12.3E4,',',ES12.3E4,')')))

    ! compute v2
    phiep(1:np) = EYE*sqrt(4.0*B1/beta(1)**2)
    phi(1:np) = cmplx(phiep,kind=DP)

    if (all(abs(phi(:)) < 2.0*(nu+1.0))) then
       call cjy(nuep,    phiep(:),J(:,1),Y(:,1))
       call cjy(nuep+1.0,phiep(:),J(:,2),Y(:,2))
    else
       do i= 1,np
          if (2.0*(nu+1.0) > abs(phi(i))) then
             call cjy(nuep,    phiep(i:i),J(i,1),Y(i,1))
             call cjy(nuep+1.0,phiep(i:i),J(i,2),Y(i,2))
          else
             call cbesj(z=phi(i),fnu=nu,kode=kode,n=num,cy=tmp(1:2),&
                  &nz=nzero,ierr=ierr)
             if (ierr > 0  .and. ierr /= 3 .and. s%quiet > 1) then
                print *, 'ERROR: CBESJ (zD=0) z=',phi(i),' nu=',nu,&
                     &' i,ierr,nz:',i,ierr,nzero
             else
                J(i,1:2) = tmp(1:2)
             end if
             
             call cbesy(z=phi(i),fnu=nu,kode=kode,n=num,cy=tmp(1:2),&
                  &nz=nzero,ierr=ierr)
             if (ierr > 0  .and. ierr /= 3 .and. s%quiet > 1) then
                print *, 'ERROR: CBESY (zD=0) z=',phi(i),' nu=',nu,&
                     &' i,ierr,nz:',i,ierr,nzero
             else
                Y(i,1:2) = tmp(1:2)
             end if
          end if
          if (any(isnan(abs(J(i,1:2)))) .or. any(isnan(abs(Y(i,1:2))))) then
             allocate(tmp2(4,0:int(nu)+2))
             call cjyva(nuep+1.0,phiep(i),nv,tmp2(1,:),tmp2(2,:),tmp2(3,:),tmp2(4,:))
             J(i,1:2) = tmp2(1,int(nv)-1:int(nv))
             Y(i,1:2) = tmp2(3,int(nv)-1:int(nv))
             deallocate(tmp2)
             print *, 'zD=0 ',nuep+1.0,phiep(i),J(i,1),Y(i,1)
          end if
       end do
    end if
    
    arg2(1:np) = beta(1)*phiep(1:np)
    aa(2,1,1:np) = arg1*J(:,1) - arg2(:)*J(:,2)
    aa(2,2,1:np) = arg1*Y(:,1) - arg2(:)*Y(:,2)

    if (s%quiet > 1) then
       write(*,998) 'zD=0  phi:',phiep(1:NPRINT),'J:',J(1:NPRINT,1),&
            &'Y:',Y(1:NPRINT,1),'a(2,1)',aa(2,1,1:NPRINT),'a(2,2)',aa(2,2,1:NPRINT)
    end if

    delta2(1:np) = aa(1,1,:)*aa(2,2,:) - aa(1,2,:)*aa(2,1,:)
    delta1(1:np) = (aa(1,1,:)*Y(:,1) - aa(1,2,:)*J(:,1))/delta2(:)* &
         & 2.0*eta(:)*sinh(eta(:)) - cosh(eta(:))

    sU(1:np,1:nz) = spread(sH(1:np,nz+1)/delta1(1:np),2,nz)*&
         &cosh(eta(1:np) .X. s%zD(1:nz))
    sD(1:np,1:nz) = sH(1:np,1:nz) + sU(1:np,1:nz)

    if (s%quiet > 1) then
       write(*,999) ' nu:',nuep,' sU:',su(1:NPRINT,1),&
            &' D1:',delta1(1:NPRINT),' D2:',delta2(1:NPRINT),' sH:',sH(1:NPRINT,1)
    end if

999 format(A,ES11.3E3,4(A,2('(',ES12.3E4,',',ES12.3E4,')')))
  end function mishraNeuman2010
  
  function mishraNeuman2010a(a,zD,s,p,f,w) result(sD)
    use constants, only : DP, EP
    use types, only : well, formation, solution
    use utility, only : operator(.X.) 
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f

    complex(EP), dimension(size(p),size(zD)) :: sD, sU
    complex(EP), dimension(size(p),size(zD)+1) :: sH
    integer ::  np, nz

    complex(EP), dimension(size(p)) :: eta, omega
    real(EP) :: lambda

    np = size(p)
    nz = size(zD)

    eta(1:np) = sqrt((a**2 + p(1:np))/f%kappa)
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0],s,p,f,w)
    lambda = f%ak - f%ac

    omega(1:np) = p(:)*exp(-f%ak*f%b1)*&
         & f%Sy*f%ac/(f%Kr*lambda)*(1.0_EP - exp(lambda*f%usL)) + &
         & a**2*f%b/f%kappa

    sU(1:np,1:nz) = spread(sH(1:np,nz+1),2,nz)*&
         & cosh(eta(1:np) .X. s%zD(1:nz))/&
         & spread(eta(:)/omega(:)*sinh(eta(:)*f%b) - &
         & cosh(eta(:)*f%b),2,nz)
    sD(1:np,1:nz) = sH(1:np,1:nz) + sU(1:np,1:nz)

  end function mishraNeuman2010a

  function mishraNeuman2010b(a,zD,s,p,f,w) result(sD)
    use constants, only : DP, EP
    use types, only : well, formation, solution
    use utility, only : operator(.X.), solve_tridiag
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f

    ! finite difference matrix related 
    ! a,b,c are sub,main,super diagonals of FD matrix
    ! sigma is result (A1 + solution in vadose zone)
    ! v is right-hand side
    complex(EP), dimension(s%order,size(p)) :: aa,b,c,sigma,v

    complex(EP), dimension(size(p),size(zD)) :: sD
    complex(EP), dimension(size(p),size(zD)+1) :: sH
    integer ::  np, nz, j, n
    integer, dimension(s%order) :: ii

    complex(EP), dimension(size(p)) :: eta, cc
    real(EP) :: lambda, h, invhsq
    complex(EP), dimension(s%order,size(p)) :: omega

    n = s%order
    h = f%usL/real(n-1,EP) ! F.D. grid spacing
    invhsq = 1.0_EP/h**2
    np = size(p)
    nz = size(zD)

    forall(j=1:n) ii(j) = j-1

    eta(1:np) = sqrt((a**2 + p(1:np))/f%kappa)
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0],s,p,f,w)
    lambda = f%ak - f%ac

    omega(1:n,1:np) = spread(p(:)*exp(-f%ak*f%b1)*&
         & f%Sy*f%ac/f%Kr,1,n)*&
         & spread(exp(lambda*(ii(1:n)*h - f%b)),2,np) + &
         & a**2/f%kappa

    ! main diagonal (first and last entries are different)
    cc(1:np) = f%ak/h - invhsq - omega(1,1:np)
    b(1,1:np) = 0.5_EP*(exp( eta(:)*f%b)*(cc(:) - eta(:)/h) + &
                      & exp(-eta(:)*f%b)*(cc(:) + eta(:)/h))
    b(2:n,1:np) = f%ak/h - 2.0*invhsq - omega(2:n,1:np)
    b(n,1:np) = b(n,1:np) + (invhsq - f%ak/h)

    ! super-diagonal (last entry (n) is undefined)
    c(1:n-1,1:np) = invhsq - f%ak/h
    c(n,1:np) = -999999.9

    ! sub-diagonal (first entry 1 is undefined, second entry is different)
    aa(2:n,1:np) = invhsq 
    aa(2,1:np) = aa(2,1:np)*cosh(eta(:)*f%b)
    aa(1,1:np) = 7777777.7

    ! right-hand side (all zero but first and second rows)
    v(3:n,1:np) = 0.0_EP
    v(2,1:np) = -sH(1:np,nz+1)*invhsq
    v(1,1:np) = (omega(1,1:np) + invhsq - f%ak/h)*sH(1:np,nz+1)

    ! sigma(1,1:np) is A1, which is constant in
    ! solution for saturated domain
    call solve_tridiag(aa,b,c,v,sigma)

    sD(1:np,1:nz) = sH(1:np,1:nz) + &
         & spread(sigma(1,:),2,nz)*&
         & cosh(eta(1:np) .X. s%zD(1:nz))

  end function mishraNeuman2010b


end module laplace_hankel_solutions




