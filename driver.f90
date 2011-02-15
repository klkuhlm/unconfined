program Driver

  ! type definitions 
  use types

  use io, only : read_input, write_timeseries_header, write_contour_header
  
  ! constants and coefficients
  use constants, only : DP, PI, EP
 
  ! function to be evaluated in Laplace/Hankel space
  use laplace_hankel_solution, only : soln => lap_hank_soln

  ! inverse Laplace transform routine (currently only de Hoog, et al)
  use inverse_Laplace_Transform, only : invlap => dehoog_invlap , pvalues => dehoog_pvalues

  implicit none

  type(invLaplace) :: l
  type(invHankel) :: h
  type(GaussLobatto) :: gl
  type(TanhSinh) :: ts
  type(well) :: w
  type(formation) :: f
  type(solution) :: s
  real(DP), parameter :: TEE_MULT = 2.0_DP
  integer, parameter :: UNIT = 20  ! output unit

  ! vectors of results and intermediate steps
  complex(EP), allocatable :: fa(:,:,:), tmp(:,:,:), GLz(:,:,:), GLarea(:,:,:)
  complex(EP), allocatable :: finint(:,:), infint(:,:), totlap(:,:)
  real(EP), allocatable :: totint(:), GLy(:)
  real(EP) :: totobs, lob, hib, width
  real(DP) :: arg
  complex(EP) :: dy
  integer :: i, j, k, m, n, nn
  integer, allocatable :: ii(:)

  ! read in data from file, do minor error checking
  ! and allocate some solution vectors
  call read_input(w,f,s,l,h,gl,ts)
  
  ! number of Laplace transform Fourier series coefficients
  l%np = 2*l%M + 1  
  ! number of tanh-sinh terms
  ts%N = 2**ts%k - 1

  allocate(finint(l%np,s%nz), infint(l%np,s%nz), totlap(l%np,s%nz), l%p(l%np), &
       & ts%w(ts%N), ts%a(ts%N), fa(ts%N,l%np,s%nz), ii(ts%N), totint(s%nz), &
       & tmp(ts%nst,l%np,s%nz), ts%kk(ts%nst), ts%NN(ts%nst), ts%hh(ts%nst))

  if (s%timeSeries) then
     call write_timeseries_header(w,f,s,l,h,gl,ts,UNIT)
  else
     call write_contour_header(w,f,s,l,h,gl,ts,UNIT)
  end if
       
  ! loop over all desired calculation times
  do i = 1, s%nt
     if (.not. s%quiet) then
        write(*,'(I4,A,ES11.4)') i,' td ',s%tD(i)
     end if
     
     ! currently using 'optimal' p values for each time
     ! TODO: compute optimal p values for a vector of times
     l%p(1:l%np) = pvalues(TEE_MULT*s%tD(i),l)

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! finite portion of Hankel integral (Tanh-Sinh quadrature)
     ! integrate from origin (a=0) to the J0 zero identified in input

     ts%kk(:) = [(ts%k-m, m=ts%nst-1,0,-1)]
     ts%NN(:) = 2**ts%kk - 1
     ts%hh(:) = 4.0_EP/(2**ts%kk)

     ! loop over all desired calculation radial distances
     do k = 1, s%nr

        ! compute abcissa and for highest-order case (others are subsets)
        ! split between finite & infinite integrals
        arg = h%j0z(h%sv(i))/s%rD(k) 
        call tanh_sinh_setup(ts,ts%k,arg)

        ! compute solution at densest set of abcissa
        !$OMP PARALLEL DO PRIVATE(nn) SHARED(fa,t,s)
        do nn = 1,ts%NN(ts%nst)
           fa(nn,1:l%np,1:s%nz) = soln(ts%a(nn),s%rD(k),l%np,s%nz,w,f,s,l)
        end do
        !$OMP END PARALLEL DO

        tmp(ts%nst,1:l%np,1:s%nz) = arg/2.0 * &
             & sum(spread(spread(ts%w(:),2,l%np),3,s%nz)*fa(:,:,:),dim=1)

        do j=ts%nst-1,1,-1
           !  only need to re-compute weights for each subsequent coarser step
           call tanh_sinh_setup(ts,ts%kk(j),arg)

           ! compute index vector, to slice up solution 
           ! for nst'th turn count regular integers
           ! for (nst-1)'th turn count even integers
           ! for (nst-2)'th turn count every 4th integer, etc...
           ii(1:ts%NN(j)) = [( m* 2**(ts%nst-j), m=1,ts%NN(j) )]

           ! estimate integral with subset of function evaluations and appropriate weights
           tmp(j,1:l%np,1:s%nz) = arg/2.0 * &
                & sum(spread(spread(ts%w(1:ts%NN(j)),2,l%np),3,s%nz) * &
                & fa(ii(1:ts%NN(j)),:,:),dim=1)
        end do
        
        if (ts%nst > 1) then
           ! perform Richardson extrapolation to spacing -> zero
           !$OMP PARALLEL DO PRIVATE(jj,k,dy) SHARED(t,tmp,finint)
           do m = 1,l%np
              do n = 1,s%nz
                 call polint(ts%hh(1:ts%nst), tmp(1:ts%nst,m,n), 0.0_EP, finint(m,n), dy)
              end do
           end do
           !$OMP END PARALLEL DO
        else
           ! no extrapolation, just use the densest (only) estimate
           finint(1:l%np,1:s%nz) = tmp(1,1:l%np,1:s%nz)
        end if

        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ! "infinite" portion of Hankel  integral for each time level
        ! integrate between zeros of Bessel function, extrapolate 
        ! area from series of areas of alternating sign
        ! TODO: use higher-order GL quadrature for lower-order zeros?

        if(.not. allocated(gl%x)) then
           allocate(gl%x(gl%ord-2), gl%w(gl%ord-2), GLy(gl%ord-2), &
                & GLarea(gl%nacc,l%np,s%nz), GLz(gl%ord-2,l%np,s%nz))
           call gauss_lobatto_setup(gl)  ! get weights and abcissa
        end if

        ! loop over regions between successive
        ! zeros of J0, integrating each independently
        do j = h%sv(i)+1, h%sv(i)+gl%nacc
           lob = real(h%j0z(j-1)/s%rD(k),EP) ! lower bound
           hib = real(h%j0z(j)/s%rD(k),EP)   ! upper bound
           width = hib - lob
           
           ! transform GL abcissa (x) to global coordinates (y)
           GLy(:) = (width*gl%x(:) + (hib+lob))/2.0
           
           !$OMP PARALLEL DO PRIVATE(k) SHARED(g,z,tD,rD,np,nz,w,f)
           do m = 1,gl%ord-2
              GLz(m,1:l%np,1:s%nz) = soln( GLy(m),s%rD(k),l%np,s%nz,w,f,s,l )
           end do
           !$OMP END PARALLEL DO
          
           GLarea(j-h%sv(i),1:l%np,1:s%nz) = width/2.0*sum(GLz(1:gl%ord-2,:,:)* &
                & spread(spread(gl%w(1:gl%ord-2),2,l%np),3,s%nz),dim=1)
        end do
 
        !$OMP PARALLEL DO PRIVATE(i,j) SHARED(area,integral)
        do j = 1,l%np
           do m = 1,s%nz
              ! accelerate each series independently
              infint(j,m) = wynn_epsilon(GLarea(1:gl%nacc,j,m))
           end do
        end do
        !$OMP END PARALLEL DO

        totlap(1:l%np,1:s%nz) = finint(1:l%np,1:s%nz) + infint(1:l%np,1:s%nz) 

        !$OMP PARALLEL DO PRIVATE(k) SHARED(totint,l,s,totlap)
        do m = 1,s%nz
           totint(m) = invlap(s%tD(i),s%tD(i)*TEE_MULT,totlap(1:l%np,m),l) 
        end do
        !$OMP END PARALLEL DO

        ! do trapezoid rule across monitoring well screen if necessary
        ! result is divided by length of interval to get average value
        if (s%timeseries) then
           if(.not. s%piezometer .and. s%zOrd > 1) then
              ! weights are uniform, except for endpoints
              ! divide by interval length
              totObs = (totint(1)/2.0 + sum(totint(2:s%zOrd)) + &
                   & totint(s%zOrd)/2.0)/(s%zOrd-1)
           else
              totObs = totint(1) ! one point, interval has no length
           end if
        end if
        
        ! write results to file
        if (s%timeSeries) then
           if (s%dimless) then
              ! dimensionless time series output
              write (UNIT,'('//s%RFMT//',1X,2('//s%HFMT//',1X))') &
                   & s%td(i), totObs ! TODO add derivative wrt logt
           else
              ! dimensional time series output
              write (UNIT,'('//s%RFMT//',1X,2('//s%HFMT//',1X))') &
                   & s%t(i), totObs*s%Hc ! TODO add derivative wrt logt
           end if
        else
           if (s%dimless) then
              ! dimensionless contour map output
              do m = 1,s%nz
                 write (UNIT,'(2('//s%RFMT//',1X),2('//s%HFMT//',1X))') &
                      & s%zD(m), s%rD(k), totint(m) ! TODO add derivative wrt logt
              end do
           else
              ! dimensional contour map output
              do m = 1,s%nz
                 write (UNIT,'(2('//s%RFMT//',1X),2('//s%HFMT//',1X))') &
                      & s%z(m), s%r(k), totint(m)*s%Hc ! TODO add derivative wrt logt
              end do
           end if
        end if
     end do
  end do

contains

  !! ###################################################
  pure subroutine tanh_sinh_setup(t,k,s)
    use constants, only : PIOV2EP
    use types, only : TanhSinh
    implicit none
    
    type(TanhSinh), intent(inout) :: t
    real(DP), intent(in) :: s
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
    P = 0.0
    xold = 2.0

    iter: do 
       if (maxval(abs(x-xold)) > spacing(1.0_EP)) then
          xold = x
          P(:,1) = 1.0
          P(:,2) = x
          
          do k = 2,N
             P(:,k+1) = ((2*k-1)*x*P(:,k) - (k-1)*P(:,k-1))/k
          end do
          
          x = xold - (x*P(:,N1)-P(:,N))/(N1*P(:,N1))
       else
          exit iter
       end if
    end do iter
    
    w = 2.0/(N*N1*P(:,N1)**2)

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
    eps(:,-1) = 0.0

    ! build up partial sums, but check for problems
    check: do i=1,np
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
    do j=0,np-2 
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

end program Driver
