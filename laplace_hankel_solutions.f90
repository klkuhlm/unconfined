!
! Copyright (c) 2012-2015 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
!

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
          ! Moench model -> just alters alphaD
          xi(1:np) = xi(:)*f%MoenchM/sum(1.0_EP/(1.0_EP + (lap%p .X. 1.0_EP/f%MoenchGamma)),dim=2)
       end if

       if (s%model == 4) then
          udp(1:np,1:nz) = theis(a,lap%p,nz)
          udp(1:np,nz+1) = udp(1:np,nz)
       else
          udp(1:np,1:nz+1) = hantush(a,[s%zD,1.0_DP],s,lap%p,f,w)
       end if

       where (spread(real(eta) < MAXEXP,2,nz))
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & cosh(eta .X. s%zD)/spread((1.0_EP + f%beta*eta*xi)*&
               & cosh(eta) + xi*sinh(eta),2,nz)
       elsewhere
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & exp(eta .X. (s%zD - 1.0_EP))/&
               & spread(1.0_EP + f%beta*eta*xi + xi,2,nz)
       end where
       deallocate(eta,xi,udp)

    case(6)
       ! Mishra/Neuman  2010 model
       select case(s%MNtype)
       case(0)
#ifdef USE_ARB_LIBRARY          
          ! naive implementation of paper (often doesnt work)
          fp(1:np,1:nz) = mishraNeuman2010(a,s%zD,s,lap%p,f,w)
#else
          stop 'ERROR: must compile code with arb dependency (and '//&
               &'USE_ARB_LIBRARY preprocessor flag) for this option'
#endif
       case(1)
          ! Malama (finiteness condition @ land surface) version 
          ! but doesn't handle partial penetration or wellbore storage
          fp(1:np,1:nz) = mishraNeumanMalama(a,s%zD,lap%p,f)
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
    real(EP) :: dD1, lD1
    integer :: np,nz

    nz = size(zD) ! # z vals requested of this function
    np = size(p)

    dD1 = 1.0_EP - w%dD
    lD1 = 1.0_EP - w%lD
    
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
     g(1,1:np,1:nz) = cosh(eta(:) .X. (dD1 - zD(1:nz)))
    ff(1,1:np,1:nz) = spread(sinh(eta(:)*w%dD),2,nz)
    ff(2,1:np,1:nz) = spread(sinh(eta(:)*lD1),2,nz)

    g(2,1:np,1:nz) = (ff(1,:,:)*cosh(eta(:) .X. zD(:)) + ff(2,:,:)*&
         & cosh(eta(:) .X. (1.0_EP - zD)))/spread(sinh(eta(:)),2,nz)

    ! below well screen
    g(3,1:np,1:nz) = spread(exp(-eta(:)*lD1) - &
         & (ff(1,:,1) + exp(-eta(:))*ff(2,:,1))/sinh(eta(:)),2,nz)

    where(zLay == 1)
       ! below well screen (0 <= zD <= 1-lD)
       udp(1:np,1:nz) = g(3,:,:)*cosh(eta(:) .X. zD(:))
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

    real(EP) :: dD1, lD1
    real(DP) :: CDw, tDb
    integer :: np,nz,i

    ! things related to Amos library
    integer(4), parameter :: kode = 1, num = 2
    integer(4) :: nzero, ierr
    integer, dimension(size(p),size(zD)) :: zLay

    np = size(p)
    nz = size(zd)

    dD1 = 1.0_EP - w%dD
    lD1 = 1.0_EP - w%lD

    ! are # zvals requested same as given in program input?
    if(sol%nz == nz) then
       zLay(1:np,1:nz) = spread(sol%zLay(:),1,np)
    else
       zLay(1:np,1:nz) = spread([sol%zLay(:),3],1,np)
    end if

    ! dimensionless pumping wellbore storage coefficient
    CDw = w%rDw**2/(2.0_DP*(w%l - w%d)*f%Ss)

    ! dimensionless well delay time
    tDb = PI*sol%rDwobs**2/(sol%sF*f%Ss)

    xi(1:np) = w%rDw*sqrt(p(:))
    eta(1:np) = sqrt((p(:) + a**2)/f%kappa)

    do i=1,np
       call cbesk(z=cmplx(xi(i),kind=DP),fnu=0.0_DP,kode=kode,&
            & n=num,cy=K(0:1,i),nz=nzero,ierr=ierr)
       if (ierr > 0 .and. ierr /= 3) then
          print *, 'ERROR: CBESK z=',xi(i),&
               &' i,ierr,nz:',i,ierr,nzero
       end if
    end do

    A0(1:np) = 2.0_EP/(p(:)*CDw*K(0,:) + xi(:)*K(1,:))
    uDf(1:np) = A0(:)/((p(:) + a**2)*(p(:)*tDb + 1.0_EP))  !!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    where (zLay == 3)
       g(1,1:np,1:nz) = cosh(eta(:) .X. (dD1 - zd(:)))
    end where

    where (zLay == 1 .or. zLay == 2)
       ff(1,1:np,1:nz) = spread(sinh(eta(:)*w%dD),2,nz)
       ff(2,1:np,1:nz) = spread(sinh(eta(:)*lD1),2,nz)
    end where

    where (zLay == 2 .or. zLay == 3)
       g(2,1:np,1:nz) = (ff(1,:,:)*cosh(eta(:) .X. zd(:)) + &
            & ff(2,:,:)*cosh(eta(:) .X. (1.0_EP - zd(:))))/&
            & spread(sinh(eta(:)),2,nz)
    end where

    where (zLay == 1)
       ff(3,1:np,1:nz) = spread(exp(-eta(:)*lD1) - &
            & (ff(1,:,1) + exp(-eta(:))*ff(2,:,1))/sinh(eta(:)),2,nz)
    end where

    where (zLay == 3) !!(spread(zD > 1-w%dD,2,nz))
       uDp(1:np,1:nz) = g(1,:,:) - g(2,:,:)
    elsewhere
       where (zLay == 1) !!(spread(zD < (1.0-w%lD),2,nz))
          uDp(1:np,1:nz) = ff(3,:,:)*cosh(eta(:) .X. zd(:))
       elsewhere
          ! zLay == 2
          uDp(1:np,1:nz) = 1.0_EP - g(2,:,:)
       end where
    end where

    u(1:np,1:nz) = spread(uDf(:)/w%bD,2,nz)*uDp(:,:)

  end function hantushstorage

#ifdef USE_ARB_LIBRARY
  function mishraNeuman2010(a,zD,s,p,f,w) result(sD)
    use constants, only : DP, EP, QP 
    use types, only : well, formation, solution
    implicit none
    
    ! https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fC_005fBINDING.html#ISO_005fC_005fBINDING
    interface
       function arb_J(nu,z) bind(c,name="arb_J") result(J)
         use, intrinsic :: iso_c_binding, only : C_FLOAT128_COMPLEX, C_FLOAT128
         complex(C_FLOAT128_COMPLEX), intent(in), value :: z
         real(C_FLOAT128), intent(in), value :: nu
         complex(C_FLOAT128_COMPLEX) :: J
       end function arb_J
    end interface
    interface
       function arb_Y(nu,z) bind(c,name="arb_Y") result(Y)
         use, intrinsic :: iso_c_binding, only : C_FLOAT128_COMPLEX, C_FLOAT128
         complex(C_FLOAT128_COMPLEX), intent(in), value :: z
         real(C_FLOAT128), intent(in), value :: nu
         complex(C_FLOAT128_COMPLEX) :: Y
       end function arb_Y
    end interface

    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f

    complex(EP), dimension(size(p),size(zD)) :: sD, sU
    complex(EP), dimension(size(p),size(zD)+1) :: sH
    integer ::  np, nz

    real(QP) :: nu, C
    complex(QP), dimension(size(p)) :: phi, arg2
    complex(QP), parameter :: EYE = (0.0_QP, 1.0_QP)

    real(QP), dimension(3) :: beta
    complex(QP) :: arg1
    complex(QP), dimension(size(p)) :: B, chi, q
    complex(QP), dimension(size(p)) :: eta
    complex(QP), dimension(size(p),0:1) :: J,Y

    integer :: i

    np = size(p)
    nz = size(zD)

    beta(1) = f%lambdaD      ! b*(ak-ac)  <lambda>
    beta(2) = f%ak*f%b1      ! ak*(psi_a - psi_k); 
    beta(3) = f%akD

    eta(1:np) = sqrt((a**2 + p(1:np))/f%kappa)
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0_DP],s,p,f,w)

    ! <between D7 & D8>
    B(1:np) = p(:)*f%ac*f%Sy/(f%Kr*f%kappa)*exp(-beta(2))
    C = (a**2)/f%kappa 

    ! <between D8 & D9>
    phi(1:np) = 2.0_QP*EYE/beta(1)*sqrt(B)*exp(0.5_QP*beta(1)*f%usLD) ! phi(L)
    nu = sqrt(beta(3)**2 + 4.0_QP*C)/beta(1) ! <n>

    do i= 1,np
       J(i,0) = arb_J(nu, phi(i))
       J(i,1) = arb_J(nu + 1.0_QP, phi(i))

       Y(i,0) = arb_Y(nu, phi(i))
       Y(i,1) = arb_Y(nu + 1.0_QP, phi(i))
    end do

    arg1 = beta(3) + nu*beta(1)
    arg2(1:np) = beta(1)*phi(1:np)

    ! <D11>
    chi(1:np) = -(arg1*J(:,0) - arg2(:)*J(:,1))/(arg1*Y(:,0) - arg2(:)*Y(:,1))

    ! phi(0)
    phi(1:np) = 2.0_QP*EYE/beta(1)*sqrt(B)

    do i= 1,np
       J(i,0) = arb_J(nu, phi(i))
       J(i,1) = arb_J(nu + 1.0_QP, phi(i))

       Y(i,0) = arb_Y(nu, phi(i))
       Y(i,1) = arb_Y(nu + 1.0_QP, phi(i))
    end do

    ! <D13>
    q(1:np) = arg1/2.0_QP - EYE*sqrt(B)*(J(:,1) + chi(:)*Y(:,1))/(J(:,0) + chi(:)*Y(:,0))

    ! <C17>  <<<<<< need to check this is properly non-dimensionalized
    sU(1:np,1:nz) = -spread(sH(1:np,nz+1)/ &
         & (cosh(eta(:)*w%bD) - eta(:)/q(:)*sinh(eta(:)*w%bD)),2,nz)* &
         & cosh(spread(eta(:),2,nz)*spread(zD(:),1,np))

    ! <B1>
    sD(1:np,1:nz) = sH(:,1:nz) + sU(:,:)

  end function mishraNeuman2010
#endif ! USE_ARB_LIBRARY
  
  function mishraNeumanMalama(a,zD,p,f) result(sD)
    use constants, only : DP, EP
    use types, only : well, formation, solution
    use utility, only : operator(.X.)
    implicit none

    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(formation), intent(in) :: f

    complex(EP), dimension(size(p),size(zD)) :: sD
    integer ::  np, nz

    real(EP) :: beta0, vartheta, phiDa, phiDk, u0
    complex(EP), dimension(size(p)) :: eta, etasq, eta1, u, v, Delta0 

    np = size(p)
    nz = size(zD)

    beta0 = f%ak*f%b  ! assume ac == ak in Malama formulation
    phiDa = f%psia/f%b
    phiDk = f%psik/f%b
    vartheta = beta0*f%Sy/(f%Ss*f%b)*exp(-beta0*(phiDa - phiDk))
    eta1(1:np) = sqrt((p(:)*vartheta + a**2)/f%kappa)

    u0 = beta0/2.0
    v(1:np) = sqrt(1.0 + (eta1(:)/u0)**2)
    u(1:np) = u0*(1.0 - v(:))

    etasq(1:np) = (p(:) + a**2)/f%kappa
    eta(1:np) = sqrt(etasq(:))
    Delta0(1:np) = eta(:)*sinh(eta(:)) - u(:)*cosh(eta(:))

    ! saturated zone solution
    sD(1:np,1:nz) = 2.0/spread(f%kappa*etasq,2,nz)* &   
         & (1.0 + spread(u/Delta0,2,nz)*cosh(eta .X. zD))

  end function mishraNeumanMalama

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
    complex(EP), dimension(size(p),s%order) :: aa,b,c,sigma,v

    complex(EP), dimension(size(p),size(zD)) :: sD
    complex(EP), dimension(size(p),size(zD)+1) :: sH
    integer ::  np, nz, j, n
    integer, dimension(s%order) :: ii

    complex(EP), dimension(size(p)) :: eta, cc
    real(EP) :: h, invhsq, B2
    real(EP), dimension(0:3) :: beta
    complex(EP), dimension(size(p),s%order) :: omega, B1
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
    sH(1:np,1:nz+1) = hantush(a,[zD,1.0_DP],s,p,f,w)

    B1(1:np,1:n) = spread(p(:)*beta(0)*exp(-beta(2))/f%kappa,2,n)
    B2 = (a**2)/f%kappa

    omega(1:np,1:n) = B1(:,:)*spread(exp(-beta(1)*ii(:)*h),1,np) + B2

    ! main diagonal (first and last entries are different)
    cc(1:np) = beta(3)/h - invhsq - omega(1:np,1)
    b(1:np,1) = 0.5_EP*(exp( eta(:))*(cc(:) - eta(:)/h) + &
                      & exp(-eta(:))*(cc(:) + eta(:)/h))
    b(1:np,2:n) = beta(3)/h - 2.0*invhsq - omega(1:np,2:n)
    b(1:np,n) = b(1:np,n) + invhsq - beta(3)/h

    ! super-diagonal (last entry (n) is undefined)
    c(1:np,1:n-1) = invhsq - beta(3)/h

    ! sub-diagonal (first entry 1 is undefined, second entry is different)
    aa(1:np,2:n) = invhsq
    aa(1:np,2) = aa(1:np,2)*cosh(eta(:))

    ! right-hand side (all zero but first and second rows)
    v(1:np,3:n) = 0.0_EP
    v(1:np,2) =   -invhsq*sH(1:np,nz+1)
    v(1:np,1) = -cc(1:np)*sH(1:np,nz+1)

    ! sigma(1:np,1) is A1, the constant needed in
    ! solution for saturated domain
    call solve_tridiag(aa,b,c,v,sigma)

    ! @ early times sigma_1 is underflowing while cosh(eta*z) is overflowing;
    ! this should allow the solution to proceed (~Hantush)
    where(spread(abs(sigma(1:np,1)) > tiny(1.0_EP),2,nz))
       sD(1:np,1:nz) = sH(1:np,1:nz) + spread(sigma(:,1),2,nz)*&
            & cosh(eta(1:np) .X. s%zD(1:nz))
    elsewhere
       sD(1:np,1:nz) = sH(1:np,1:nz)
    end where

    if (s%quiet > 1) then
       write(*,999) ' sD:',sD(1:NPRINT,1), ' A1:',sigma(1:NPRINT,1),&
            & ' a01:',aa(1:NPRINT,1),' b01:',b(1:NPRINT,1),' c01:',c(1:NPRINT,1)&
            & ,' v01:',v(1:NPRINT,1),' omega01:',omega(1:NPRINT,1)
       do j=2,n
          write(*,998) &
               &'                                                   '//&
               &'                                                          a',j,':'&
               &,aa(1:NPRINT,j),' b',j,':',b(1:NPRINT,j),' c',j,':',c(1:NPRINT,j)&
               & ,' v',j,':',v(1:NPRINT,j),' omega',j,':',omega(1:NPRINT,j)
       end do
    end if

998 format(5(A,I2.2,A,2('(',ES11.2E4,',',ES11.2E4,')')))
999 format(7(A,2('(',ES11.2E4,',',ES11.2E4,')')))

  end function mishraNeuman2010FD

end module laplace_hankel_solutions

