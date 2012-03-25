module utility
  implicit none
  private
  public :: logspace, linspace, is_finite, operator(.X.), solve_tridiag, spec_basis

  interface operator(.X.)
     module procedure outerprod_zd, outerprod_dz, outerprod_dd, outerprod_zz
  end interface

contains

  function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    integer :: i
    real(DP) :: rnum, range, sgn

    if (num == 1) then
       ! convention puts a single point in the
       ! middle of the interval
       v = [(lo + hi)/2.0]
    else
       rnum = real(num - 1,DP)
       range = abs(hi - lo)
       sgn = sign(1.0_DP,hi-lo) ! if lo > high, count backwards
       forall (i=0:num-1)
          v(i+1) = lo + sgn*real(i,DP)*range/rnum
       end forall
    end if
  end function linspace

  function logspace(lo,hi,num) result(v)
    use constants, only : DP
    integer, intent(in) :: lo,hi,num
    real(DP), dimension(num) :: v
    v = 10.0_DP**linspace(real(lo,DP),real(hi,DP),num)
  end function logspace

  elemental function is_finite(x) result(pred)
    use constants, only : EP
    complex(EP), intent(in) :: x
    logical :: pred
    pred = .not. (isnan(abs(x)) .or. abs(x) > huge(abs(x)))
  end function is_finite

  pure function outerprod_dd(da,db) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da,db
    real(DP), dimension(size(da),size(db)) :: c
    c = spread(da,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(da))
  end function outerprod_dd

  pure function outerprod_zd(za,db) result(c)
    use constants, only : EP,DP
    complex(EP), intent(in), dimension(:) :: za
    real(DP), intent(in), dimension(:) :: db
    complex(EP), dimension(size(za),size(db)) :: c
    c = spread(za,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(za))
  end function outerprod_zd

  pure function outerprod_dz(da,zb) result(c)
    use constants, only : EP,DP
    real(DP), intent(in), dimension(:) :: da
    complex(EP), intent(in), dimension(:) :: zb
    complex(EP), dimension(size(da),size(zb)) :: c
    c = spread(da,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(da))
  end function outerprod_dz

  pure function outerprod_zz(za,zb) result(c)
    use constants, only : EP
    complex(EP), intent(in), dimension(:) :: za,zb
    complex(EP), dimension(size(za),size(zb)) :: c
    c = spread(za,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(za))
  end function outerprod_zz

  pure subroutine solve_tridiag(a,b,c,v,x)
    use constants, only : EP
    implicit none
    ! modified from Wikipedia
    ! http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm

    !      a - sub-diagonal (means it is the diagonal below the main diagonal)
    !      b - the main diagonal
    !      c - sup-diagonal (means it is the diagonal above the main diagonal)
    !      v - right part
    !      x - the answer
    !      n - number of equations

    complex(EP),dimension(:,:),intent(in) :: a,b,c,v
    complex(EP),dimension(size(a,1),size(a,2)),intent(out) :: x
    complex(EP),dimension(size(a,1),size(a,2)) :: bp,vp
    complex(EP),dimension(size(a,2)) :: m
    integer :: i,np,n

    n = size(a,1)
    np = size(a,2)

    ! Make copies of the b and v variables so that they are unaltered by this sub
    bp(1,1:np) = b(1,:)
    vp(1,1:np) = v(1,:)

    !The first pass (setting coefficients):
    firstpass: do i = 2,n
       m(1:np) = a(i,:)/bp(i-1,:)
       bp(i,1:np) = b(i,:) - m(:)*c(i-1,:)
       vp(i,1:np) = v(i,:) - m(:)*vp(i-1,:)
    end do firstpass

    x(n,1:np) = vp(n,:)/bp(n,:)
    !The second pass (back-substition)
    backsub:do i = n-1, 1, -1
       x(i,1:np) = (vp(i,:) - c(i,:)*x(i+1,:))/bp(i,:)
    end do backsub

  end subroutine solve_tridiag

  subroutine spec_basis(x,nbasis,phi,phix,phixx)
    use constants, only : EP

    integer, intent(in) :: nbasis
    real(EP), intent(in) :: x
    real(EP), intent(out), dimension(nbasis) :: phi,phix,phixx

    integer :: i,n,nn
    real(EP) :: t,c,s,tnt,tntt

    intrinsic :: sign

    ! compute the basis functions for spectral problem
    ! adapted from Boyd, "Chebyshev and Fourier Spectral
    ! Methods", 2000, section 6.15 and Appendix A.2

    if (abs(x) < 1.0) then
       t = acos(x)
       c = cos(t)
       s = sin(t)
       do i=1,nbasis
          n = i+1
          phi(i) = cos(n*t)
          tnt = -n*sin(n*t)
          tntt = -n*n*phi(i)

          ! convert t-derivatives into x-derivatives
          phix(i) = -tnt/s  ! x-derivative of nth Cheb poly
          phixx(i) = tntt/(s*s) - tnt*c/(s*s*s) ! second x-derivative
       end do
    else
       ! alternative formulation when |x|=1
       do i=1,nbasis
          n = i+1
          nn = n**2
          phi(i) = sign(1.0_EP,x)**n
          phix(i) =  sign(1.0_EP,x)**(n+1) * nn
          phixx(i) = sign(1.0_EP,x)**n * nn * (nn-1.0)/3.0
       end do
    end if

  end subroutine spec_basis

end module utility


