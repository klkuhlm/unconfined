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

    complex(EP), allocatable :: eta(:), xi(:), udp(:,:)
    integer :: mm

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

    case(3:5)
      ! Moench hybrid water table condition made a special case
      ! Malama 2011 partial penetrating model (Neuman 74 when beta=0)
       allocate(eta(np),xi(np),udp(np,nz+1))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)

       if (s%model == 3) then
          mm = f%MoenchAlphaM
          ! Moench model -> just alters alphaD
          xi(1:np) = xi(:)*mm/sum(1.0/(1.0 + (lap%p .X. 1.0/f%MoenchGamma)),dim=2)
       end if

       if (s%model == 4) then
          udp(1:np,1:nz) = theis(a,lap%p,nz)
          udp(1:np,nz+1) = udp(1:np,nz)
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
       select case(s%MNtype)
       case(0)
          ! naive implementation of paper (often doesnt work)
          fp(1:np,1:nz) = mishraNeuman2010(a,s%zD,s,lap%p,f,w)
       case(2)
          ! finite difference solution of ODE in vadose zone
          fp(1:np,1:nz) = mishraNeuman2010FD(a,s%zD,s,lap%p,f,w)
       end select

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
    use constants, only : DP, EP
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
    use constants, only : DP, EP, PI
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

    ! dimensionless pumping wellbore storage coefficient
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
    use constants, only : DP, EP, EYE
    use types, only : well, formation, solution
    use utility, only : operator(.X.)
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

    real(EP) :: nuep, B2 !, nv
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

    do i= 1,np
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
    end do

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

    do i= 1,np
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
    end do

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

  function mishraNeuman2010FD(a,zD,s,p,f,w) result(sD)
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
    real(EP) :: h, invhsq, B2
    real(EP), dimension(0:3) :: beta
    complex(EP), dimension(s%order,size(p)) :: omega, B1
    integer, parameter :: NPRINT = 2

    n = s%order
    h = f%usLD/real(n-1,EP) ! F.D. grid spacing
    invhsq = 1.0_EP/h**2
    np = size(p)
    nz = size(zD)

    ! integers 0:n-1
    forall(j=1:n) ii(j) = j-1

    beta(0) = f%ac*f%Sy/f%Ss
    beta(1) = -f%lambdaD  ! b*(ac-ak)  (-1* definition in driver_io.f90)
    beta(2) = f%ak*f%b1   ! ak*(psi_a - psi_k)
    beta(3) = f%akD       ! ak*b

    eta(1:np) = sqrt((a**2 + p(1:np))/f%kappa)
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0],s,p,f,w)

    B1(1:n,1:np) = spread(p(:)*beta(0)*exp(-beta(2))/f%kappa,1,n)
    B2 = (a**2)/f%kappa

    omega(1:n,1:np) = B1(:,:)*spread(exp(-beta(1)*ii(:)*h),2,np) + B2

    ! main diagonal (first and last entries are different)
    cc(1:np) = beta(3)/h - invhsq - omega(1,1:np)
    b(1,1:np) = 0.5_EP*(exp( eta(:))*(cc(:) - eta(:)/h) + &
                      & exp(-eta(:))*(cc(:) + eta(:)/h))
    b(2:n,1:np) = beta(3)/h - 2.0*invhsq - omega(2:n,1:np)
    b(n,1:np) = b(n,1:np) + invhsq - beta(3)/h

    ! super-diagonal (last entry (n) is undefined)
    c(1:n-1,1:np) = invhsq - beta(3)/h
!!$    c(n,1:np) = -999999.9

    ! sub-diagonal (first entry 1 is undefined, second entry is different)
    aa(2:n,1:np) = invhsq
    aa(2,1:np) = aa(2,1:np)*cosh(eta(:))
!!$    aa(1,1:np) = 7777777.7

    ! right-hand side (all zero but first and second rows)
    v(3:n,1:np) = 0.0_EP
    v(2,1:np) =   -invhsq*sH(1:np,nz+1)
    v(1,1:np) = -cc(1:np)*sH(1:np,nz+1)

    ! sigma(1,1:np) is A1, the constant needed in
    ! solution for saturated domain
    call solve_tridiag(aa,b,c,v,sigma)

    ! @ early times sigma_1 is underflowing while cosh(eta*z) is overflowing;
    ! this should allow the solution to proceed (~Hantush)
    where(spread(abs(sigma(1,1:np)) > tiny(1.0_EP),2,nz))
       sD(1:np,1:nz) = sH(1:np,1:nz) + spread(sigma(1,:),2,nz)*&
            & cosh(eta(1:np) .X. s%zD(1:nz))
    elsewhere
       sD(1:np,1:nz) = sH(1:np,1:nz)
    end where

    if (s%quiet > 1) then
       write(*,999) ' sD:',sD(1:NPRINT,1), ' A1:',sigma(1,1:NPRINT),&
            & ' a01:',aa(1,1:NPRINT),' b01:',b(1,1:NPRINT),' c01:',c(1,1:NPRINT)&
            & ,' v01:',v(1,1:NPRINT),' omega01:',omega(1,1:NPRINT)
       do j=2,n
          write(*,998) &
               &'                                                   '//&
               &'                                                          a',j,':'&
               &,aa(j,1:NPRINT),' b',j,':',b(j,1:NPRINT),' c',j,':',c(j,1:NPRINT)&
               & ,' v',j,':',v(j,1:NPRINT),' omega',j,':',omega(j,1:NPRINT)
       end do
    end if

998 format(5(A,I2.2,A,2('(',ES11.2E4,',',ES11.2E4,')')))
999 format(7(A,2('(',ES11.2E4,',',ES11.2E4,')')))

  end function mishraNeuman2010FD

end module laplace_hankel_solutions

