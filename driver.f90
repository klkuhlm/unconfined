Program Driver

! needed for intel compiler
#ifdef INTEL
  use ifport, only : dbesj0, dbesj1     
#endif

  ! Pass data to routines using modules
  use shared_data         ! utility.f90

  ! constants and coefficients
  use constants, only : DP, PI, EP
 
  use lap_hank_soln, only : unconfined_wellbore_slug
  use inverse_Laplace_Transform, only : dehoog_invlap, dehoog_pvalues 

  use utilities, only : logspace

  implicit none

  real(DP), allocatable :: j0zero(:)
  real(DP) :: x, dx

  character(75) :: outfile, infile, timefname
  integer :: i, numt, terms, minlogsplit, maxlogsplit, np
  integer :: splitrange, zerorange
  real(EP), allocatable :: GLx(:),GLw(:)
  integer :: minlogt,maxlogt
  
  ! vectors of results and intermediate steps
  integer, allocatable :: splitv(:)
  real(DP), allocatable  :: ts(:), td(:), dt(:) 
  complex(EP), allocatable :: fa(:,:)
  complex(EP), allocatable :: finint(:), infint(:), totlap(:), p(:)
  real(EP) :: totint
  real(DP) :: tee

  !! k parameter in tanh-sinh quadrature (order is 2^k - 1)
  integer :: tsk, tsN, nst
  complex(EP), allocatable :: tmp(:,:)
  real(EP), allocatable :: tsw(:), tsa(:), tshh(:)
  integer, allocatable :: tskk(:), tsNN(:), ii(:)
  integer :: GLorder, j0split(1:2), naccel, err, m, nn, jj, j, numtfile, numtcomp

  ! dimensional things
  logical :: quiet, dimless, computetimes
  real(DP) :: B, Kr, sigma, Ss, Sy
  real(DP) :: l,d,LL,LLe,rw,rc,g,nu,Tc
  complex(EP) :: polerr

  intrinsic :: get_command_argument

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! read input parameters from file some minimal error checking
  call get_command_argument(1,infile)
  if (len_trim(infile) == 0) then
!!$     write(*,'(A)') 'no command-line input filename supplied, using default input filename: input.dat'
     infile = 'input.dat'
  end if

  open(unit=19, file=infile, status='old', action='read',iostat=err)
  if(err /= 0) then
     write(*,'(A)') 'ERROR opening main input file '//trim(infile)//' for reading'
     stop
  end if
  
  read(19,*) quiet, dimless          ! suppress output to screen, write dimensionless output?
  read(19,*) B              ! initial saturated aquifer thickness [L]
  read(19,*) l, d           ! distance from aquifer top to bottom & top of packer (l > d) [L]
  read(19,*) rw, rc         ! well and casing radii [L]
  read(19,*) Kr, kappa         ! aquifer radial K [L^2/T] and anisotropy ratio (Kz/Kr)
  read(19,*) Ss, Sy         ! aquifer specific storage [1/L] and specific yield [-]
  read(19,*) g, nu          ! grav accel [L/T^2], kinematic viscosity [L^2/T]

  if (.not. quiet) then
     write(*,'(A,2(L1,1X))') 'quiet?, dimless output? ', quiet, dimless
     write(*,'(A,ES14.7)') 'B (initial aquier sat thickness):',B  
     write(*,'(A,2(ES14.7,1X))') 'l,d (packer bot & top from top of aquifer):',l,d  
     write(*,'(A,2(ES14.7,1X))') 'rw,rc (well and casing radii):',rw,rc  
     write(*,'(A,2(ES14.7,1X))') 'Kr,kappa',Kr, kappa
     write(*,'(A,2(ES14.7,1X))') 'Ss,Sy',Ss, Sy
     write(*,'(A,2(ES14.7,1X))') 'g,nu (grav acc & kin viscosity in consistent units)',g,nu 
  end if

  if(any([B,Kr,kappa,Ss,l,rw,rc,g,nu] < epsilon(1.0))) then
     write(*,'(A)') 'input error: zero or negative distance or aquifer parameters' 
     stop
  end if

  if (d >= l) then
     write(*,'(2(A,ES11.5))') 'input error: l must be > d; l=',l,' d=',d
     stop
  end if
  
  read(19,*) lap%M, lap%alpha, lap%tol

  if (lap%M < 1) then
     write(*,'(A,I0)')  'deHoog # FS terms must be > 0 M=',lap%M
     stop 
  end if

  if (lap%tol < epsilon(lap%tol)) then
     lap%tol = epsilon(lap%tol)
     write(*,'(A,ES11.5)') 'WARNING: increased INVLAP solution tolerance to ',lap%tol 
  end if

  read(19,*) tsk, nst                       ! tanh-sinh quadrature parameters
  read(19,*) j0split(1:2), naccel, GLorder  ! Gauss-Lobatto quadrature parameters

  if(tsk - nst < 2) then
     write(*,'(2(A,I0),A)') 'Tanh-Sinh k is too low (',tsk,') for given level&
          & of Richardson extrapolation (',nst,').  Increase tsk or decrease nst.'
     stop
  end if  

  if(any([j0split(:),naccel, tsk] < 1)) then
     write(*,'(A)') 'max/min split, # accelerated terms and tanh-sinh k must be >= 1' 
     stop
  end if

  terms = maxval(j0split(:))+naccel+1
  allocate(j0zero(terms))

  read(19,*) minlogt,maxlogt,numtcomp           ! log_10(t_min), log_10(t_max), # times
  read(19,*) computetimes,numtfile,timefname    ! compute times?, # times in file, time file  
  read(19,*) outfile                            ! output filename
  
  if (computetimes) then
     numt = numtcomp
  else
     numt = numtfile
  end if
  
  ! solution vectors
  allocate(finint(2*lap%M+1), infint(2*lap%M+1), totlap(2*lap%M+1), p(2*lap%M+1),&
       & ts(numt), td(numt), splitv(numt), dt(numt))

  if (computetimes) then
     ! computing times 
     ts = logspace(minlogt,maxlogt,numtcomp)
  else
     open(unit=22, file=timefname, status='old', action='read',iostat=err)
     if(err /= 0) then
        write(*,'(A)') 'ERROR opening time data input file '//trim(timefname)//' for reading'
        stop 
     else
        ! times listed one per line
        do i=1,numtfile
           read(22,*,iostat=err) ts(i)
           if(err /= 0) then
              write(*,'(A,I0,A)') 'ERROR reading time ',i,' from input file '//trim(timefname)
              stop
           end if
        end do
     end if
     close(22)
  end if
  
  ! compute derived or dimensionless properties
  gamma = 1.0   ! no skin, always 1
  sigma = Sy/(Ss*B)
  alphaD = kappa/sigma
  rDw = rw/B
  bD = (l-d)/B
  lD = l/B
  dD = d/B
  Tc = B**2/(Kr/Ss)
  LL =    d + (l-d)/2.0*(rc/rw)**4
  LLe =  LL + (l-d)/2.0*(rc/rw)**2
  beta(1) = 8.0*nu*LL/(rc**2 *g*Tc)
  beta(2) = LLe/(g*Tc**2)
  CD = (rc/B)**2/(Ss*(l-d))
  td(1:numt) = ts(:)/Tc

  if (.not. quiet) then
     write(*,'(A,ES14.7)') 'kappa:   ',kappa
     write(*,'(A,ES14.7)') 'sigma:   ',sigma
     write(*,'(A,ES14.7)') 'alpha_D: ',alphaD
     write(*,'(A,ES14.7)') 'r_{D,w}: ',rDw
     write(*,'(A,ES14.7)') 'b_D: ',bD
     write(*,'(A,ES14.7)') 'l_D: ',lD
     write(*,'(A,ES14.7)') 'Tc:  ',Tc
     write(*,'(A,ES14.7)') 'L:   ',LL
     write(*,'(A,ES14.7)') 'Le:  ',LLe
     write(*,'(A,2(ES14.7,1X))') 'beta: ',beta
     write(*,'(A,ES14.7)') 'C_D: ',CD
     write(*,'(A,I0,2(ES14.7,1X))') 'deHoog: M,alpha,tol: ',lap%M, lap%alpha, lap%tol
     write(*,'(A,2(I0,1X))'), 'tanh-sinh: k, num extrapolation steps ',tsk,nst
     write(*,'(A,4(I0,1X))'), 'GL quad: J0 split, num zeros accel, GL-order ',j0split(:),naccel,GLorder
     write(*,'(A,L1)') 'compute times? ',computetimes
     if(computetimes) then
        write(*,'(A,3(I0,1X))') 'log10(tmin), log10(tmax), num times ',minlogt,maxlogt,numtcomp
     else
        write(*,'(A,I0,1X,A)') 'num times, filename for t inputs ',numtfile,trim(timefname)
     end if
  end if

  ! compute zeros of J0 bessel function
  !************************************
  ! asymptotic estimate of zeros - initial guess
  j0zero = [((real(i-1,DP) + 0.75_DP)*PI,i=1,terms)]
  !$OMP PARALLEL DO PRIVATE(i,x,dx) SHARED(j0zero,TERMS)
  do i=1,TERMS
     x = j0zero(i)
     NR: do
        ! Newton-Raphson
        dx = dbesj0(x)/dbesj1(x)
        x = x + dx
        if(abs(dx) < spacing(x)) then
           exit NR
        end if
     end do NR
     j0zero(i) = x
  end do
  !$OMP END PARALLEL DO

  ! split between finite/infinite part should be 
  ! small for large time, large for small time
  do i=1,numt
     zerorange = maxval(j0split(:)) - minval(j0split(:))
     minlogsplit = floor(minval(log10(td)))   ! min log(td) -> maximum split
     maxlogsplit = ceiling(maxval(log10(td))) ! max log(td) -> minimum split
     splitrange = maxlogsplit - minlogsplit + 1 
     splitv(i) = minval(j0split(:)) + &
          & int(zerorange*((maxlogsplit - log10(td(i)))/splitrange))
  end do
 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! finite portion of Hankel integral for each time level
  
  np = 2*lap%M+1  ! number of Laplace transform Fourier series coefficients

  if(.not. allocated(tsw)) then
     tsN = 2**tsk - 1
     allocate(tsw(tsN),tsa(tsN),fa(tsN,np),ii(tsN),tmp(nst,np))
     allocate(tskk(nst),tsNN(nst),tshh(nst))
  end if

  ! open output file or die
  open (unit=20, file=outfile, status='replace', action='write', iostat=err)
  if (err /= 0) then
     print *, 'cannot open output file for writing' 
     stop
  end if
  
  ! echo input parameters at head of output file
  write(20,'(A,L1)') '# quiet? ', quiet
  write(20,'(A,ES14.7)') '# B (initial aquier sat thickness): ',B  
  write(20,'(A,2(ES14.7,1X))') '# l,d (packer bot & top from top of aquifer): ',l,d  
  write(20,'(A,2(ES14.7,1X))') '# rw,rc (well and casing radii): ',rw,rc  
  write(20,'(A,2(ES14.7,1X))') '# Kr,kappa: ',Kr, kappa
  write(20,'(A,2(ES14.7,1X))') '# Ss,Sy: ',Ss, Sy
  write(20,'(A,2(ES14.7,1X))') '# g,nu (grav acc & kin viscosity in consistent units): ',g,nu 
  write(20,'(A,I0,2(ES14.7,1X))') '# deHoog M, alpha, tol: ',lap%M,lap%alpha,lap%tol
  write(20,'(A,2(I0,1X))'), '# tanh-sinh: k, num extrapolation steps: ',tsk,nst
  write(20,'(A,4(I0,1X))'), '# GL quad: J0 split, # zeros accel, GL-order: ',j0split(:),naccel,GLorder
  write(20,'(A,L1)') '# compute times?: ',computetimes
  if(computetimes) then
     write(20,'(A,3(I0,1X))') '# log10(tmin), log10(tmax), num times: ',minlogt,maxlogt,numtcomp
  else
     write(20,'(A,I0,1X,A)') '# num times, filename for t inputs: ',numtfile,trim(timefname)
  end if
  if (dimless) then
     write (20,'(A,/,A,/,A)') '#','#         t_D                       H^{-1}[ L^{-1}[ f_D ]] ',&
          & '#------------------------------------------------------------'
  else
     write (20,'(A,/,A,/,A)') '#','#         t                         H^{-1}[ L^{-1}[ f ]]   ',&
          & '#------------------------------------------------------------'
  end if
     
  ! loop over all requested times
  do i = 1, numt
     if (.not. quiet) then
        write(*,'(I4,A,ES11.4)') i,' td ',td(i)
     end if
     tsval = td(i)

     ! slug-test system
     !===========================

     tskk(1:nst) = [(tsk-m, m=nst-1,0,-1)]
     tsNN(1:nst) = 2**tskk - 1
     tshh(1:nst) = 4.0_EP/(2**tskk)

     ! compute abcissa and for highest-order case (others are subsets)
     call tanh_sinh_setup(tskk(nst),j0zero(splitv(i))/rDw,tsw(1:tsNN(nst)),tsa(1:tsNN(nst)))

     ! compute solution at densest set of abcissa
     !$OMP PARALLEL DO PRIVATE(nn) SHARED(fa,tsa,tsNN,nst)
     do nn=1,tsNN(nst)
        fa(nn,1:np) = unconfined_wellbore_slug(tsa(nn))
     end do
     !$OMP END PARALLEL DO

     !$OMP PARALLEL WORKSHARE
     tmp(nst,1:np) = (j0zero(splitv(i))/rDw)/2.0_DP*sum(spread(tsw(:),2,np)*fa(:,:),dim=1)
     !$OMP END PARALLEL WORKSHARE

     do j=nst-1,1,-1
        !  only need to re-compute weights for each subsequent step
        call tanh_sinh_setup(tskk(j),j0zero(splitv(i))/rDw,tsw(1:tsNN(j)),tsa(1:tsNN(j)))

        ! compute index vector, to slice up solution 
        ! for nst'th turn count regular integers
        ! for (nst-1)'th turn count even integers
        ! for (nst-2)'th turn count every 4th integer, etc...
        ii(1:tsNN(j)) = [( m* 2**(nst-j), m=1,tsNN(j) )]

        ! estimate integral with subset of function evaluations and appropriate weights
        tmp(j,1:np) = (j0zero(splitv(i))/rDw)/2.0_DP* &
             & sum(spread(tsw(1:tsNN(j)),2,np)*fa(ii(1:tsNN(j)),1:np),dim=1)
     end do

     if (nst > 1) then
        !$OMP PARALLEL DO PRIVATE(jj,PolErr) SHARED(tshh,tmp,finint)
        do jj=1,np
           call polint(tshh(1:nst),tmp(1:nst,jj),0.0_EP,finint(jj),PolErr)
        end do
        !$OMP END PARALLEL DO
     else
        finint(:) = tmp(1,1:np)
     end if


     !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! "infinite" portion of Hankel  integral for each time level
     ! integrate between zeros of Bessel function, extrapolate 
     ! area from series of areas of alternating sign

     if(.not. allocated(GLx)) then
        allocate(Glx(GLorder-2),GLw(GLorder-2))
        call gauss_lobatto_setup(GLorder,GLx,GLw)  ! get weights and abcissa
     end if

     tsval = td(i)

     infint(1:np) = inf_integral(unconfined_wellbore_slug,splitv(i),naccel,&
          &glorder,j0zero,GLx,GLw)

     ! perform numerical inverse Laplace transform and sum infinite and finite
     ! portions of Hankel integral
     tee = td(i)*2.0
     p = dehoog_pvalues(tee,lap)
     totlap(1:np) = finint(1:np) + infint(1:np) ! omega
     totlap(1:np) = beta(1) + beta(2)*p + gamma*totlap(1:np)/2.0 ! f
     totlap(1:np) = totlap(:)/(1.0 + p(:)*totlap(:))  ! Psi
     totint = dehoog_invlap(td(i),tee,totlap(1:np),lap) 

     ! write results to file (as we go)
     if (dimless) then
        write (20,'(3(1x,ES24.15E3))') td(i), totint
     else
        write (20,'(3(1x,ES24.15E3))') ts(i), totint
     end if
     
  end do
  deallocate(finint,infint,ts,td,splitv,dt)

contains

  !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
  !!   end of program flow
  !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !! ###################################################
  function inf_integral(integrand,split,num,ord,zeros,GLx,GLw) result(integral)
    use constants, only : DP,EP
    use shared_data, only : rdw, lap
    implicit none

    interface 
       function integrand(a) result(H)
         use constants, only : EP
         use shared_data, only : lap
         real(EP), intent(in) :: a
         complex(EP), dimension(2*lap%M+1) :: H
       end function integrand
    end interface

    integer, intent(in) :: split
    integer, intent(in) :: num,ord
    real(DP), dimension(:), intent(in) :: zeros

    real(EP), dimension(ord-2) :: GLy
    complex(EP), dimension(ord-2,2*lap%M+1) :: GLz
    real(EP), dimension(:), intent(in) :: GLx,GLw
    complex(EP), dimension(num,2*lap%M+1) :: area
    real(EP) :: lob,hib,width
    complex(EP), dimension(2*lap%M+1) :: integral
    integer :: i, j, k, np

    np = 2*lap%M+1

    do j = split+1, split+num
       lob = real(zeros(j-1)/rdw,EP) ! lower bound
       hib = real(zeros(j)/rdw,EP)   ! upper bound
       width = real(hib,EP) - real(lob,EP)
       
       ! transform GL abcissa to global coordinates
       GLy(:) = (width*GLx(:) + (hib + lob))/2.0
          
       !$OMP PARALLEL DO PRIVATE(k) SHARED(GLz,GLy,ord,np)
       do k=1,ord-2
          GLz(k,1:np) = integrand( GLy(k) )
       end do
       !$OMP END PARALLEL DO
          
       area(j-split,1:np) = width/2.0* sum(GLz(1:ord-2,:)*spread(GLw(1:ord-2),2,np),dim=1)
    end do
 
    !$OMP PARALLEL DO PRIVATE(i) SHARED(area,integral,num)
    do i=1,np
       ! accelerate each series independently
       integral(i) = wynn_epsilon(area(1:num,i))
    end do
    !$OMP END PARALLEL DO
    
  end function inf_integral

  !! ###################################################
  subroutine tanh_sinh_setup(k,s,w,a)
    use constants, only : PIOV2EP
    implicit none
    
    integer, intent(in) :: k
    real(DP), intent(in) :: s
    real(EP), intent(out), dimension(2**k-1) :: w, a

    integer :: N,r,i
    real(EP) :: h
    real(EP), allocatable :: u(:,:)

    !! compute weights 
    N = 2**k-1
    r = (N-1)/2
    h = 4.0_EP/2**k
    allocate(u(2,N))
       
    u(1,1:N) = PIOV2EP*cosh(h*[(i, i=-r,r)])
    u(2,1:N) = PIOV2EP*sinh(h*[(i, i=-r,r)])
    
    a(1:N) = tanh(u(2,1:N))
    w(1:N) = u(1,1:N)/cosh(u(2,:))**2
    w(1:N) = 2.0_EP*w(1:N)/sum(w(1:N))

    ! map the -1<=x<=1 interval onto 0<=a<=s
    a = (a + 1.0_EP)*s/2.0_EP

  end subroutine tanh_sinh_setup


  !! ###################################################
  subroutine gauss_lobatto_setup(ord,abcissa,weight)
    use constants, only : PIEP, EP
    implicit none
    integer, intent(in) :: ord
    real(EP), intent(out), dimension(ord-2) :: weight, abcissa
    real(EP), dimension(ord,ord) :: P
    real(EP), dimension(ord) :: x,xold,w
    integer :: i, N, N1, k

    ! leave out the endpoints (abcissa = +1 & -1), since 
    ! they will be the zeros of the Bessel functions
    ! (therefore, e.g., 5th order integration only uses 3 points)
    ! code copied from Matlab routine by Greg von Winckel, at
    ! http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4775

    N = ord-1
    N1 = N+1
    ! first guess
    x = cos(PIEP*[(i, i=0,N)]/N)
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

    ! leave off endpoints
    abcissa = x(2:ord-1)
    weight = w(2:ord-1)

  end subroutine gauss_lobatto_setup


  !! ###################################################
  !! wynn-epsilon acceleration of partial sums, given a series
  !! all intermediate sums / differences are done in extended precision
  function wynn_epsilon(series) result(accsum)
    use constants, only : EONE, EP, DP

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
       if((abs(series(i)) > huge(0.0D0)).or.(series(i) /= series(i))) then
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
          if(abs(denom) > tiny(EONE)) then ! check for div by zero
             eps(m,j+1) = eps(m+1,j-1) + EONE/denom
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
  SUBROUTINE polint(xa,ya,x,y,dy)
    ! xa and ya are given x and y locations to fit an nth degree polynomial
    ! through.  results is a value y at given location x, with error estimate dy

    use constants, only : DP, EP
    IMPLICIT NONE
    REAL(EP), DIMENSION(:), INTENT(IN) :: xa
    complex(EP), dimension(:), intent(in) :: ya
    REAL(EP), INTENT(IN) :: x
    complex(EP), INTENT(OUT) :: y,dy
    INTEGER :: m,n,ns
    COMPLEX(EP), DIMENSION(size(xa)) :: c,d,den
    REAL(EP), DIMENSION(size(xa)) :: ho
    integer, dimension(1) :: tmp

    n=size(xa)
    c=ya
    d=ya
    ho=xa-x
    tmp = minloc(abs(x-xa))
    ns = tmp(1)
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(abs(den(1:n-m)) == 0.0)) then
          print *, 'polint: calculation failure'
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
  END SUBROUTINE polint

end program Driver
