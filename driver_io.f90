module io
implicit none

private
public :: read_input, write_header

contains
  subroutine read_input(w,f,s,lap,h,gl,ts)
    use constants, only : DP, PI
    use types
    use utilities, only : logspace
    
    ! only needed for intel compiler
#ifdef INTEL
    use ifport, only : dbesj0, dbesj1     
#endif

    type(invLaplace), intent(inout) :: lap
    type(invHankel), intent(inout) :: h
    type(GaussLobatto), intent(inout) :: gl
    type(TanhSinh), intent(inout) :: ts
    type(well), intent(inout) :: w
    type(formation), intent(inout) :: f
    type(solution), intent(inout) :: s

    integer :: ioerr, terms
    character(39) :: fmt

    ! intermediate steps in vsplit
    integer :: zerorange, minlogsplit, maxlogsplit, splitrange

    ! used to compute times
    integer :: numtfile, numtcomp, minlogt, maxlogt

    character(NUMCHAR) :: timeFileName, spaceFileName
    real(DP) :: tval, rval

    ! used to compute J0 zeros
    real(DP) :: x, dx 

    intrinsic :: get_command_argument

    ! read input parameters from file some minimal error checking
    call get_command_argument(1,infile)
    if (len_trim(infile) == 0) then
       ! default input file name if no command-line argument
       infile = 'input.dat'
    end if

    open(unit=19, file=infile, status='old', action='read',iostat=ioerr)
    if(ioerr /= 0) then
       write(*,'(A)') 'ERROR opening main input file '//trim(infilename)//&
            & ' for reading'
       stop
    end if
  
    ! ## switches and settings #####

    ! suppress output to screen,  which model to use?, 
    ! write dimensionless output?, 
    ! time-series (compute one location through time)? 
    !     [non-time-series option is contours at a single time],
    ! piezometer (point monitoring well) or integrate finite observation screen length?,
    read(19,*) s%quiet, s%model, s%dimless, s%timeseries, s%piezometer

    ! ## pumping well parameters #####

    ! volumetric pumping rate [L^3/T]
    read(19,*) w%Q

    ! distance from aquifer top to bottom & top of packer (l > d) [L]
    read(19,*) w%l, w%d

    ! pumping well: well and casing radii [L]
    read(19,*) w%rw, w%rc 

    ! dimensionless skin
    read(19,*) f%gamma

    ! pumping well time behavior
    read(19,'(I3)', advance='no') lap%timeType 
    if (lap%timeType > -1) then
       ! functional time behavior (two or fewer parameters)
       allocate(lap%timePar(2))
       read(19,*) lap%timePar(:)
       if (.not. s%quiet) then
          write(*,'(A)') lap%timeDescrip(lap%timeType)
          write(*,'(A,2(1X,'//s%rfmt//'))') 'time behavior, par1, par2 ::', &
               & lap%timeType, lap%timePar(:)
       end if
    else
       ! arbitrarily long piecewise-constant time behavior 
       allocate(lap%timePar(-2*lap%timeType+1))
       read(19,*) lap%timePar(:)
       if (.not. s%quiet) then
          fmt = '(A,    ('//s%rfmt//',1X),A,    ('//s%rfmt//',1X))'
          write(fmt(8:11),'(I4.4)')  size(lap%timePar(:-lap%timeType+1),1)
          write(fmt(26:29),'(I4.4)') size(lap%timePar(-lap%timeType+2:),1)
          write(*,'(A)') lap%timeDescrip(9)
          write(*,fmt) 'time behavior:  ti, tf | Q multiplier each step :: ', &
               & lap%timeType,lap%timePar(:-lap%timeType+1),'| ',&
               & lap%timePar(-lap%timeType+2:)
       end if
    end if

    ! ## aquifer / formation parameters #####

    ! initial saturated aquifer thickness [L]
    read(19,*) f%b

    ! aquifer radial K [L^2/T] and anisotropy ratio (Kz/Kr)
    read(19,*) f%Kr, f%kappa
    
    ! aquifer specific storage [1/L] and specific yield [-]
    read(19,*) f%Ss, f%Sy

    ! unsaturated zone thickness [L] and exponential sorbtive parameter [1/L]
    read(19,*) f%usl, f%usalpha 
    
    ! ## echo check parameters #####

    if (.not. s%quiet) then
       write(*,'(A,L1,1X,A,1X)') 'model, dimless output?:',&
            & s%quiet, s%modelDescription(s%model) ,s%dimless
       write(*,'(A,2(L1,1X))') 'hydrograph?, piezometer?:', &
            & s%hydrograph, s%piezometer
       write(*,'(A,'//s%rfmt//')') 'b (initial aquier sat thickness):',f%b  
       write(*,'(A,2('//s%rfmt//',1X))') 'l,d (screen bot&top from above):', w%l, w%d  
       write(*,'(A,2('//s%rfmt//',1X))') 'rw,rc (well and casing radii):', w%rw, w%rc  
       write(*,'(A,3('//s%rfmt//',1X))') 'Kr,kappa,gamma:', f%Kr, f%kappa, f%gamma
       write(*,'(A,2('//s%rfmt//',1X))') 'Ss,Sy:', f%Ss, f%Sy
    end if

    if(any([f%Sy, f%gamma] < 0.0)) then
       write(*,*) 'ERROR: negative aquifer parameters:', &
            &  f%Sy, f%gamma
       stop
    end if   

    ! spacing(0.) is the spacing between computer representable numbers @ 0
    if(any([w%b,f%Kr,f%kappa,f%Ss,w%l,w%rw,w%rc,us%L] < spacing(0.0))) then
       write(*,*) 'ERROR: zero or negative parameters:', &
            &  w%b, f%Kr, f%kappa, f%Ss, w%l, w%rw, w%rc
       stop
    end if

    if (w%d >= w%l) then
       write(*,'(2(A,'//s%rfmt//'))') 'ERROR: l must be > d; l=',w%l,' d=',w%d
       stop
    end if
  
    ! ## numerical implementation-related parameters #####

    ! Laplace transform (deHoog et al) parameters
    read(19,*) lap%M, lap%alpha, lap%tol

    ! tanh-sinh quadrature parameters
    ! integration order 2^(k-1) and Richardson extrapolation order
    read(19,*) ts%k, ts%nst
    
    ! Gauss-Lobatto quadrature parameters
    ! max/min J0 zero to split at, # zeros to accelerate, GL-order
    read(19,*) h%j0s(1:2), gl%nacc, gl%ord  

    ! ## checking of numerical parameters #####

    if (lap%M < 2) then
       write(*,'(A,I0)')  'ERROR: deHoog # FS terms must be >= 1 M=',lap%M
       stop 
    end if

    if (lap%tol < epsilon(lap%tol)) then
       lap%tol = epsilon(lap%tol)
       write(*,'(A,'//s%rfmt//')') 'WARNING: increased INVLAP tolerance to ',&
            & lap%tol 
    end if

    if(ts%k - ts%nst < 2) then
       write(*,'(2(A,I0),A)') 'ERROR: Tanh-Sinh k is too low (',ts%k,&
            & ') for given level of Richardson extrapolation (',ts%nst,&
            &').  Increase k or decrease nst.'
       stop
    end if

    if(any([h%j0split(:),gl%naccel, ts%k] < 1)) then
       write(*,'(2A,4(I0,1X))') 'ERROR max/min split, # accelerated terms, ',&
            & 'and tanh-sinh k must be >= 1:',h%j0split(:),gl%naccel, ts%k
       stop
    end if

    ! ## locations and times where solution is desired #####

    ! solution at one location through time
    ! times can be specified in a file,
    ! or log-sampled between endpoints

    read(19,*) timeFileName, tval  ! either tval or data in file used
    
    ! solution at one or more locations through time
    ! locations can be specified in a file,
    ! or linear-sampled between endpoints

    read(19,*) spaceFileName, rval ! either rval or data in file used

    ! if computing contour or profile, these aren't used
    if (s%piezometer) then
       ! point observation depth
       read(19,*) w%pztop
    else
       ! observation location with a finite screen
       read(19,*) w%pztop, w%pzbot
    end if

    if (s%timeseries) then
       ! if computing a time series, read time-related
       ! parameters from input file

       allocate(w%r(1)) ! only one observation distance
       w%r = rval

       open(unit=22, file=trim(timeFileName), status='old', &
            & action='read',iostat=ioerr)
       if(ioerr /= 0) then
          write(*,'(A)')) 'ERROR opening time input file '// &
               & trim(timeFileName)//' for reading'
       end if
       
       ! compute times?, # times to read (not used if vector computed)
       read(22,*) computeTimes, numtfile

       ! log_10(t_min), log_10(t_max), # times
       read(22,*) minlogt, maxlogt, numtcomp   

       if (computeTimes) then
          s%nt = numtcomp
       else
          s%nt = numtfile
       end if

       allocate(s%t(s%nt), s%tD(s%nt), s%sv(s%nt))

       if (s%computetimes) then
          ! computing times 
          s%t = logspace(minlogt,maxlogt,numtcomp)
       else
          ! times listed one per line
          do i=1,numtfile
             read(22,*,iostat=ioerr) s%t(i)
             if(ioerr /= 0) then
                write(*,'(A,I0,A)') 'ERROR reading time ',i,&
                     & ' from input file '//trim(timefilename)
                stop
             end if
          end do
          close(22)
       end if
       
    else
       ! if computing a contour map or profile
       ! read space related parameters from other file

       allocate()

    end if
    
    ! output filename
    read(19,*) s%outfilename                        
    close(19)

    ! characteristic length
    s%Lc = f%b

    ! characteristic time 
    s%Tc = f%b**2/(f%Kr/f%Ss)
    
    ! compute derived or dimensionless properties
    f%sigma = f%Sy/(f%Ss*f%b)
    f%alphaD = f%kappa/f%sigma

    ! dimensionless lengths
    w%rDw = w%rw/s%Lc
    w%bD = (w%l - w%d)/s%Lc
    w%lD = w%l/s%Lc
    w%dD = w%d/s%Lc

    ! dimensionless times
    s%tD(:) = s%t(:)/s%Tc

    if (.not. s%quiet) then
       write(*,'(A,'//s%rfmt//')') 'kappa:   ',f%kappa
       write(*,'(A,'//s%rfmt//')') 'sigma:   ',f%sigma
       write(*,'(A,'//s%rfmt//')') 'alpha_D: ',f%alphaD
       write(*,'(A,'//s%rfmt//')') 'Tc:  ',s%Tc
       write(*,'(A,'//s%rfmt//')') 'Lc:  ',s%Lc
       write(*,'(A,'//s%rfmt//')') 'r_{D,w}: ',w%rDw
       write(*,'(A,'//s%rfmt//')') 'b_D: ',w%bD
       write(*,'(A,'//s%rfmt//')') 'l_D: ',w%lD
       write(*,'(A,I0,2('//s%rfmt//',1X))') 'deHoog: M,alpha,tol: ', &
            & lap%M, lap%alpha, lap%tol
       write(*,'(A,2(I0,1X))'), 'tanh-sinh: k, num extrapolation steps ', &
            & ts%k, ts%nst
       write(*,'(A,4(I0,1X))'), 'GL: J0 split, num zeros accel, GL-order ',&
            & h%j0split(:), gl%naccel, gl%order
       write(*,'(A,L1)') 'compute times? ',s%computetimes
       if(s%computetimes) then
          write(*,'(A,3(I0,1X))') 'log10(tmin), log10(tmax), num times ',&
               & minlogt, maxlogt, numtcomp
       else
          write(*,'(A,I0,1X,A)') 'num times, filename for t inputs ',&
               & numtfile, trim(timefilename)
       end if
    end if

    terms = maxval(h%j0s(:)) + gl%nacc + 1
    allocate(h%j0z(terms))

    ! compute zeros of J0 bessel function
    ! (quicker / more accurate than reading from file)
    !************************************
    ! asymptotic estimate of zeros - initial guess
    forall (i=0:terms-1) h%j0zero(i+1) = (i + 0.75)*PI
    do i=1,terms
       x = h%j0zero(i)
       NR: do
          ! Newton-Raphson
          dx = dbesj0(x)/dbesj1(x)
          x = x + dx
          if(abs(dx) < spacing(x)) then
             exit NR
          end if
       end do NR
       h%j0zero(i) = x
    end do
    
    allocate(splitv(s%numt))

    ! split between finite/infinite part should be 
    ! small for large time, large for small time
    zerorange = maxval(h%j0split(:)) - minval(h%j0split(:))
    minlogsplit = floor(minval(log10(s%tD)))   ! min log(td) -> maximum split
    maxlogsplit = ceiling(maxval(log10(s%tD))) ! max log(td) -> minimum split
    splitrange = maxlogsplit - minlogsplit + 1 
    splitv = minval(h%j0split(:)) + &
         & int(zerorange*((maxlogsplit - log10(s%tD))/splitrange))

  end subroutine read_input

  subroutine write_header(w,f,s,lap,h,gl,ts,unit)
    use types
    
    type(invLaplace), intent(in) :: lap
    type(invHankel), intent(in) :: h
    type(GaussLobatto), intent(in) :: gl
    type(TanhSinh), intent(in) :: ts
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    type(solution), intent(in) :: s
    integer, intent(in) :: unit
    
    integer :: ioerr

    ! open output file or die
    open (unit=unit, file=s%outfilename, status='replace', action='write', iostat=ioerr)
    if (ioerr /= 0) then
       write(*,'(A)') 'cannot open output file '//trim(s%outfilename)//' for writing' 
       stop
    end if
  
    ! echo input parameters at head of output file
    write(20,'(A,L1)') '# quiet?, dimensionless? ::', &
         & s%quiet, s%dimless
    write(20,'(A,'//s%rfmt//')') '# b (initial aquier sat thickness) ::', &
         & f%b
    write(20,'(A,2('//s%rfmt//',1X))') '# l,d (screen bot & top) ::',&
         & w%l, w%d  
    write(20,'(A,2('//s%rfmt//',1X))') '# rw,rc (well/casing radii) ::',&
         & w%rw, w%rc  
    write(20,'(A,2('//s%rfmt//',1X))') '# Kr,kappa ::', f%Kr, f%kappa
    write(20,'(A,2('//s%rfmt//',1X))') '# Ss,Sy :: ',f%Ss, f%Sy
    write(20,'(A,'//s%rfmt//')') '# gamma :: ',f%gamma
    write(20,'(A,I0,2('//s%rfmt//',1X))') '# deHoog M, alpha, tol ::',&
         & lap%M, lap%alpha, lap%tol
    write(20,'(A,2(I0,1X))'), '# tanh-sinh: k, n extrapolation steps ::',&
         & ts%k, ts%nst
    write(20,'(A,4(I0,1X))'), '# GLquad: J0 split, n 0-accel, GL-order ::',&
         & h%j0split(:), gl%naccel, gl%order
    write(20,'(A,L1)') '# times ::',s%numt
    if (s%dimless) then
       write (20,'(A,/,A,/,A)') '#','#         t_D                       '//&
            & 'H^{-1}[ L^{-1}[ f_D ]] ',&
            & '#------------------------------------------------------------'
    else
       write (20,'(A,/,A,/,A)') '#','#         t                         '//&
            & 'H^{-1}[ L^{-1}[ f ]]   ',&
            & '#------------------------------------------------------------'
    end if

  end subroutine write_header
end module io

