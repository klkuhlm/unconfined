module integration
implicit none

private
public :: tanh_sinh_setup, gauss_lobatto_setup, wynn_epsilon, polint

contains
  !! ###################################################
  subroutine tanh_sinh_setup(t,k,s,j)
    use constants, only : PIOV2EP, DP, EP
    use types, only : TanhSinh
    implicit none
    
    type(TanhSinh), intent(inout) :: t
    real(EP), intent(in) :: s
    integer, intent(in) :: k, j

    integer :: N,r,i
    real(EP) :: h
    real(EP), allocatable :: u(:,:)

    !! compute weights 
    N = 2**k-1
    r = (N-1)/2
    h = 4.0_EP/2**k
    allocate(u(2,N))
     
    forall (i=-r:r)
       u(1,i+r+1) = PIOV2EP*cosh(h*i)
       u(2,i+r+1) = PIOV2EP*sinh(h*i)
    end forall

    t%Q(j)%w(1:N) = u(1,:)/cosh(u(2,:))**2
    t%Q(j)%w(1:N) = 2.0_EP*t%Q(j)%w(:)/sum(t%Q(j)%w(:))

    ! only compute abcissa, if vector is allocated
    if (allocated(t%Q(j)%a)) then
       ! TODO: only use half interval with bunched abcissa @ origin?
       ! map the -1<=x<=1 interval onto 0<=a<=s
       t%Q(j)%a(1:N) = (tanh(u(2,:)) + 1.0_EP)*s/2.0_EP
    end if

    deallocate(u)

  end subroutine tanh_sinh_setup

  !! ###################################################
  subroutine gauss_lobatto_setup(gl)
    use constants, only : PIEP, EP
    use types, only : GaussLobatto
    implicit none

    type(GaussLobatto), intent(inout) :: gl
    real(EP), dimension(gl%ord,gl%ord) :: P
    real(EP), dimension(gl%ord) :: x, xold, w
    integer :: i, N, N1, k

    ! leave out the endpoints (abcissa = +1 & -1), since 
    ! they will be the zeros of the Bessel functions
    ! (therefore, e.g., 5th order integration only uses 3 points)
    ! code modified from Matlab routine by Greg von Winckel, at
    ! http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4775

    N = gl%ord-1
    N1 = N+1

    ! first guess
    forall (i=0:N) 
       x(i+1) = cos(PIEP*i/N)
    end forall

    ! initialize Vandermonde matrix
    P = 0.0_EP
    xold = 2.0_EP

    iter: do 
       if (maxval(abs(x-xold)) > spacing(1.0_EP)) then
          xold = x
          P(:,1) = 1.0_EP
          P(:,2) = x
          
          do k = 2,N
             P(:,k+1) = ((2*k-1)*x*P(:,k) - (k-1)*P(:,k-1))/k
          end do
          
          x = xold - (x*P(:,N1)-P(:,N))/(N1*P(:,N1))
       else
          exit iter
       end if
    end do iter
    
    w = 2.0_EP/(N*N1*P(:,N1)**2)

    ! leave off endpoints (BF defined as zero there)
    gl%x = x(2:gl%ord-1)
    gl%w = w(2:gl%ord-1)

  end subroutine gauss_lobatto_setup

  !! ###################################################
  !! wynn-epsilon acceleration of partial sums, given a series
  !! all intermediate sums / differences are done in extended precision
  function wynn_epsilon(series,quiet) result(accsum)
    use constants, only : EP
    use utility, only : is_finite
    implicit none

    integer, parameter :: MINTERMS = 4
    complex(EP), dimension(:), intent(in) :: series
    integer, intent(in) :: quiet
    complex(EP) :: accsum, denom 
    integer :: ns, i, j, m
    complex(EP), dimension(1:size(series),-1:size(series)-1) :: eps

    ns = size(series)

    ! build up partial sums, but check for problems
    check: do i=1,ns
!!$       write(*,'(I0,A,ES12.3E4,A,ES12.3E4,A)',advance='no'), i,'(',real(series(i)),',',aimag(series(i)),') '
       if (.not. is_finite(series(i))) then
          ns = i-1
          if(ns < MINTERMS) then
             if (quiet > 1) then
                write(*,'(A)',advance='no') 'not enough Wynn-Epsilon series to accelerate '
             end if
             accsum = -999999.9  ! make it clear answer is bogus
             goto 777
          else
             if (quiet > 1) then
                write(*,'(A,I3,A)',advance='no') 'Wynn-Epsilon series&
                     &, truncated to ',ns,' terms. '
             end if
             exit check
          end if
       else
          ! term is good, continue
          eps(i,0) = sum(series(1:i))
       end if
    end do check

    ! first column is intiallized to zero
    eps(:,-1) = 0.0_EP

    ! build up epsilon table (each column has one less entry)
    do j = 0,ns-2 
       do m = 1,ns-(j+1)
          denom = eps(m+1,j) - eps(m,j)
          if(abs(denom) > epsilon(abs(denom))) then ! check for div by zero
             eps(m,j+1) = eps(m+1,j-1) + 1.0_EP/denom
          else
             accsum = eps(m+1,j)
             if (quiet > 1) then
                write(*,'(A,I0,1X,I0,A)') 'epsilon cancel ',m,j,':'
!!$                write(*,'(A,I0,1X,I0,3(A,ES12.3E4,A,ES12.3E4),A)',advance='no') &
!!$                     & 'epsilon cancel ',m,j,' denom(',real(denom),',',aimag(denom),&
!!$                     &') eps(m+1,j)=(',real(eps(m+1,j)),',',aimag(eps(m+1,j)),&
!!$                     &') eps(m,j)=(',real(eps(m,j)),',',aimag(eps(m,j)),') : '
             end if
             goto 777
          end if
       end do
    end do

    ! if made all the way through table use "corner value" of triangle as answer
    if(mod(ns,2) == 0) then
       accsum = eps(2,ns-2)  ! even number of terms - corner is acclerated value
    else
       accsum = eps(2,ns-3)  ! odd numbers, use one column in from corner
    end if
777 continue

  end function wynn_epsilon

  !! ###################################################
  ! polynomial extrapolation modified from numerical recipes f90 (section 3.1)
  subroutine polint(xa,ya,x,y,dy)
    ! xa and ya are given x and y locations to fit an nth degree polynomial
    ! through.  results is a value y at given location x, with error estimate dy

    ! x is real and extended-precision
    ! y is complex and extended-precision

    use constants, only : DP, EP
    implicit none

    real(EP), dimension(:), intent(IN) :: xa
    complex(EP), dimension(:), intent(in) :: ya
    real(EP), intent(IN) :: x
    complex(EP), intent(OUT) :: y,dy
    integer :: m,n,ns
    complex(EP), dimension(size(xa)) :: c,d,den
    real(EP), dimension(size(xa)) :: ho

    n=size(xa)
    c=ya
    d=ya
    ho=xa-x
    ns = sum(minloc(abs(x-xa))) ! sum turns 1-element vector to a scalar
    y=ya(ns)
    ns=ns-1
    ! TODO vectorize this loop
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(abs(den(1:n-m)) < spacing(0.0))) then
          write(*,*) 'polint: calculation failure',abs(den(1:n-m))
          stop
       end if
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  end subroutine polint
end module integration
