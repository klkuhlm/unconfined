module integration
implicit none

private
public :: tanh_sinh_setup, gauss_lobatto_setup, wynn_epsilon, polint

contains
  !! ###################################################
  subroutine tanh_sinh_setup(t,k,s)
    use constants, only : PIOV2EP, DP, EP
    use types, only : TanhSinh
    implicit none
    
    type(TanhSinh), intent(inout) :: t
    real(EP), intent(in) :: s
    integer, intent(in) :: k

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
    
    if(size(t%a) /= N) then
       deallocate(t%a,t%w)
       allocate(t%a(N),t%w(N))
    end if
    
    t%a(1:N) = tanh(u(2,1:N))
    t%w(1:N) = u(1,:)/cosh(u(2,:))**2

    deallocate(u)
    t%w(1:N) = 2.0_EP*t%w(:)/sum(t%w(:))

    ! TODO: only use half interval with bunched abcissa @ origin?
    ! map the -1<=x<=1 interval onto 0<=a<=s
    t%a = (t%a + 1.0_EP)*s/2.0_EP

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
    forall (i=0:N) x(i+1) = cos(PIEP*i/N)
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
  function wynn_epsilon(series) result(accsum)
    use constants, only : EP, DP

    implicit none
    integer, parameter :: MINTERMS = 4
    complex(EP), dimension(1:), intent(in) :: series
    complex(DP) :: accsum
    complex(EP) :: denom
    integer :: np, i, j, m
    complex(EP), dimension(1:size(series),-1:size(series)-1) :: eps

    np = size(series)

    ! first column is intiallized to zero
    eps(:,-1) = 0.0_EP

    ! build up partial sums, but check for problems
    check: do i = 1,np
       if((abs(series(i)) > huge(1.0_EP)).or.(series(i) /= series(i))) then
          ! +/- Infinity or NaN ? truncate series to avoid corrupting accelerated sum
          np = i-1
          if(np < MINTERMS) then
             write(*,'(A)',advance='no') 'not enough Wynn-Epsilon series to accelerate '
             accsum = -999999.9  ! make it clear answer is bogus
             goto 777
          else
             write(*,'(A,I3,A)',advance='no') 'Wynn-Epsilon series&
                  &, truncated to ',np,' terms. '
             exit check
          end if
       else
          ! term is good, continue
          eps(i,0) = sum(series(1:i))
       end if
    end do check
    
    ! build up epsilon table (each column has one less entry)
    do j = 0,np-2 
       do m = 1,np-(j+1)
          denom = eps(m+1,j) - eps(m,j)
          if(abs(denom) > spacing(0.0)) then ! check for div by zero
             eps(m,j+1) = eps(m+1,j-1) + 1.0/denom
          else
             accsum = real(eps(m+1,j),DP)
             write(*,'(A,I2,1X,I2,A)',advance='no') 'epsilon cancel ',m,j,' '
             goto 777
          end if
       end do
    end do

    ! if made all the way through table use "corner value" of triangle as answer
    if(mod(np,2) == 0) then
       accsum = real(eps(2,np-2),DP)  ! even number of terms - corner is acclerated value
    else
       accsum = real(eps(2,np-3),DP)  ! odd numbers, use one column in from corner
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
