program Driver

  ! type definitions 
  use types

  use driver_io, only : read_input, write_timeseries_header, write_contour_header
  
  ! constants and coefficients
  use constants, only : DP, PI, EP, RFMT, HFMT
 
  ! function to be evaluated in Laplace/Hankel space
  use laplace_hankel_solutions, only : soln => lap_hank_soln

  ! inverse Laplace transform routine (currently only de Hoog, et al)
  use invLap, only : dehoog => dehoog_invlap , pvalues => dehoog_pvalues

  ! integration setup and extrapolation routines
  use integration, only : tanh_sinh_setup, gauss_lobatto_setup, wynn_epsilon, polint

  ! openMP library interface 
  !$ use omp_lib, only : OMP_get_max_threads, OMP_get_num_procs

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
  real(EP) :: totobs, lob, hib, width, arg
  complex(EP) :: dy
  integer :: i, j, k, m, n, nn
  integer, allocatable :: iv(:)

  !$ write(*,'(2(A,I0))') '# avail. processors: ',OMP_get_num_procs(), &
  !$     & '  maximum # threads: ',OMP_get_max_threads()

  ! read in data from file, do minor error checking
  ! and allocate some solution vectors
  call read_input(w,f,s,l,h,gl,ts)
  
  ! number of Laplace transform Fourier series coefficients
  l%np = 2*l%M + 1  
  ! number of tanh-sinh terms
  ts%N = 2**ts%k - 1

  allocate(finint(l%np,s%nz), infint(l%np,s%nz), totlap(l%np,s%nz), l%p(l%np), &
       & ts%w(ts%N), ts%a(ts%N), fa(ts%N,l%np,s%nz), iv(ts%N), totint(s%nz), &
       & tmp(ts%Rord,l%np,s%nz), ts%kv(ts%Rord), ts%Nv(ts%Rord), ts%hv(ts%Rord))

  if (s%timeSeries) then
     call write_timeseries_header(w,f,s,l,h,gl,ts,UNIT)
  else
     call write_contour_header(w,f,s,l,h,gl,ts,UNIT)
  end if
       
  ! loop over all desired calculation times
  do i = 1, s%nt
     if (.not. s%quiet) then
        write(*,'(I5,A,'//RFMT//')') i,' td ',s%tD(i)
     end if
     
     ! using 'optimal' vector of p values for each time
     l%p(1:l%np) = pvalues(TEE_MULT*s%tD(i),l)
!!$     write(*,'(A,7(A,ES9.3,A,ES9.3,A))') 'p(1:7) ',&
!!$          &('(',real(l%p(k)),',',aimag(l%p(k)),')',k=1,7)

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! finite portion of Hankel integral (Tanh-Sinh quadrature)
     ! integrate from origin (a=0) to the J0 zero identified in input

     forall (m = 1:ts%Rord) 
        ! count up to k, a vector of length Rord
        ts%kv(m) = ts%k - ts%Rord + m
     end forall
     ts%Nv(:) = 2**ts%kv - 1
     ts%hv(:) = 4.0_EP/(2**ts%kv)

     ! loop over all desired calculation radial distances
     do k = 1, s%nr

        ! compute abcissa and for highest-order case (others are subsets)
        ! split between finite & infinite integrals
        arg = h%j0z(h%sv(i))/s%rD(k) 
        call tanh_sinh_setup(ts,ts%k,arg)

        ! compute solution at densest set of abcissa
        !$OMP PARALLEL DO SHARED(fa,ts,w,f,s,l,k)
        do nn = 1,ts%Nv(ts%Rord)
           fa(nn,1:l%np,1:s%nz) = soln(ts%a(nn),s%rD(k),l%np,s%nz,w,f,s,l)
        end do
        !$OMP END PARALLEL DO

        tmp(ts%Rord,1:l%np,1:s%nz) = arg/2.0 * &
             & sum(spread(spread(ts%w(:),2,l%np),3,s%nz)*fa(:,:,:),dim=1)

        do j=ts%Rord-1,1,-1
           !  only need to re-compute weights for each subsequent coarser step
           call tanh_sinh_setup(ts,ts%kv(j),arg)

           ! compute index vector, to slice up solution 
           ! for nst'th turn count regular integers
           ! for (nst-1)'th turn count even integers
           ! for (nst-2)'th turn count every 4th integer, etc...
           iv(1:ts%Nv(j)) = [( m* 2**(ts%Rord-j), m=1,ts%Nv(j) )]

           ! estimate integral with subset of function evaluations and appropriate weights
           tmp(j,1:l%np,1:s%nz) = arg/2.0 * &
                & sum(spread(spread(ts%w(1:ts%Nv(j)),2,l%np),3,s%nz) * &
                & fa(iv(1:ts%Nv(j)),:,:),dim=1)
        end do
        
        if (ts%Rord > 1) then
           ! perform Richardson extrapolation to spacing -> zero
           !$OMP PARALLEL DO SHARED(ts,tmp,finint)
           do m = 1,l%np
              do n = 1,s%nz
                 call polint(ts%hv(1:ts%Rord), tmp(1:ts%Rord,m,n), 0.0_EP, finint(m,n), dy)
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
           
           !$OMP PARALLEL DO SHARED(GLy,w,f,s,l,k)
           do m = 1,gl%ord-2
              GLz(m,1:l%np,1:s%nz) = soln( GLy(m),s%rD(k),l%np,s%nz,w,f,s,l )
           end do
           !$OMP END PARALLEL DO
          
           GLarea(j-h%sv(i),1:l%np,1:s%nz) = width/2.0*sum(GLz(1:gl%ord-2,:,:)* &
                & spread(spread(gl%w(1:gl%ord-2),2,l%np),3,s%nz),dim=1)
        end do
 
        !$OMP PARALLEL DO SHARED(GLarea,infint)
        do j = 1,l%np
           do m = 1,s%nz
              ! accelerate each series independently
              infint(j,m) = wynn_epsilon(GLarea(1:gl%nacc,j,m))
           end do
        end do
        !$OMP END PARALLEL DO

        totlap(1:l%np,1:s%nz) = finint(1:l%np,1:s%nz) + infint(1:l%np,1:s%nz) 

        !$OMP PARALLEL DO SHARED(totint,l,s,totlap,i)
        do m = 1,s%nz
           totint(m) = dehoog(s%tD(i),s%tD(i)*TEE_MULT,totlap(1:l%np,m),l) 
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
              write (UNIT,'('//RFMT//',1X,2('//HFMT//',1X))') &
                   & s%td(i), totObs ! TODO add derivative wrt logt
           else
              ! dimensional time series output
              write (UNIT,'('//RFMT//',1X,2('//HFMT//',1X))') &
                   & s%t(i), totObs*s%Hc ! TODO add derivative wrt logt
           end if
        else
           if (s%dimless) then
              ! dimensionless contour map output
              do m = 1,s%nz
                 write (UNIT,'(2('//RFMT//',1X),2('//HFMT//',1X))') &
                      & s%zD(m), s%rD(k), totint(m) ! TODO add derivative wrt logt
              end do
           else
              ! dimensional contour map output
              do m = 1,s%nz
                 write (UNIT,'(2('//RFMT//',1X),2('//HFMT//',1X))') &
                      & s%z(m), s%r(k), totint(m)*s%Hc ! TODO add derivative wrt logt
              end do
           end if
        end if
     end do
  end do

end program Driver
