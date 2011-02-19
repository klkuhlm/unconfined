module driver_io
implicit none

private
public :: read_input, write_timeseries_header, write_contour_header

contains
  subroutine read_input(w,f,s,lap,h,gl,ts)
    use constants, only : EP, DP, PI, PIEP, NUMCHAR, RFMT, HFMT
    use types, only : invLaplace, invHankel, GaussLobatto, tanhSinh, well, formation, solution
    use utility, only : logspace, linspace
    
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
    integer :: zrange, minlsp, maxlsp, spRange

    ! used to compute times/spaces
    integer :: numtfile, numtcomp, minlogt, maxlogt
    integer :: numRComp, numZComp, numRFile, numZFile
    real(DP) :: minR, maxR, minZ, maxZ

    character(NUMCHAR) :: timeFileName, spaceFileName, inputFileName
    logical :: computeSpace, computeTimes
    real(DP) :: tval, rval
    integer :: i

    ! used to compute J0 zeros
    real(EP) :: x, dx 

    intrinsic :: get_command_argument, bessel_j0, bessel_j1

    ! read input parameters from file some minimal error checking
    call get_command_argument(1,inputFileName)
    if (len_trim(inputFileName) == 0) then
       ! default input file name if no command-line argument
       inputFileName = 'input.dat'
    end if

    open(unit=19, file=inputFileName, status='old', action='read',iostat=ioerr)
    if(ioerr /= 0) then
       write(*,'(A)') 'ERROR opening main input file '//trim(inputFileName)//&
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
    else
       ! arbitrarily long piecewise-constant time behavior 
       allocate(lap%timePar(-2*lap%timeType+1))
       read(19,*) lap%timePar(:)
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
       write(*,'(A,A,1X,1L)') 'model, dimless output?:: ',&
            & trim(s%modelDescrip(s%model)) ,s%dimless
       write(*,'(A,2(L1,1X))') 'time-series plot?, piezometer?:: ', &
            & s%timeSeries, s%piezometer
       write(*,'(A,'//RFMT//')') 'Q:: ',w%Q
       if (lap%timeType > -1) then
          write(*,'(A)') 'pumping well time behavior :: '//trim(lap%timeDescrip(lap%timeType))
          write(*,'(A,I0,2(1X,'//RFMT//'))') 'time behavior, par1, par2 ::', &
               & lap%timeType, lap%timePar(:)
       else
          fmt = '(A,    ('//RFMT//',1X),A,    ('//RFMT//',1X))'
          write(fmt(8:11),'(I4.4)')  size(lap%timePar(:-lap%timeType+1),1)
          write(fmt(26:29),'(I4.4)') size(lap%timePar(-lap%timeType+2:),1)
          write(*,'(A)') 'pumping well time behavior :: '//trim(lap%timeDescrip(9))
          write(*,fmt) 'time behavior:  ti, tf | Q multiplier each step :: ', &
               & lap%timeType,lap%timePar(:-lap%timeType+1),'| ',&
               & lap%timePar(-lap%timeType+2:)
       end if
       write(*,'(A,'//RFMT//')') 'b (initial aquier sat thickness):: ',f%b  
       write(*,'(A,2('//RFMT//',1X))') 'l/b,d/b (relative screen bot&top from above):: ', w%l, w%d  
       write(*,'(A,2('//RFMT//',1X))') 'rw,rc (well and casing radii):: ', w%rw, w%rc  
       write(*,'(A,3('//RFMT//',1X))') 'Kr,kappa,gamma:: ', f%Kr, f%kappa, f%gamma
       write(*,'(A,2('//RFMT//',1X))') 'Ss,Sy:: ', f%Ss, f%Sy
    end if

    if(any([f%Sy, f%gamma, w%l] < 0.0)) then
       write(*,*) 'ERROR: negative aquifer parameters:', &
            &  f%Sy, f%gamma
       stop
    end if   

    ! spacing(0.) is the spacing between computer representable numbers @ 0
    if(any([f%b,f%Kr,f%kappa,f%Ss,w%rw,w%rc] < spacing(0.0))) then
       write(*,*) 'ERROR: zero or negative parameters:', &
            &  f%b, f%Kr, f%kappa, f%Ss, w%l, w%rw, w%rc
       stop
    end if

    if (w%d <= w%l) then
       write(*,'(2(A,'//RFMT//'))') 'ERROR: d must be > l; l=',w%l,' d=',w%d
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
       write(*,'(A,'//RFMT//')') 'WARNING: increased INVLAP tolerance to ',&
            & lap%tol 
    end if

    if(ts%k - ts%nst < 2) then
       write(*,'(2(A,I0),A)') 'ERROR: Tanh-Sinh k is too low (',ts%k,&
            & ') for given level of Richardson extrapolation (',ts%nst,&
            &').  Increase k or decrease nst.'
       stop
    end if

    if(any([h%j0s(:),gl%nacc, ts%k] < 1)) then
       write(*,'(2A,4(I0,1X))') 'ERROR max/min split, # accelerated terms, ',&
            & 'and tanh-sinh k must be >= 1:',h%j0s(:),gl%nacc, ts%k
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
    ! point observation depth
    ! only top used if piezometer
    read(19,*) s%zTop, s%zBot, s%zOrd 

    if (s%zTop >= s%zBot .or. s%zTop < 0.0 .or. s%zBot > 1.0) then
       write(*,*) 'ERROR: top of monitoring well screen must be',&
            & 'above bottom and both between 0 and 1',s%zTop,s%zBot
       stop
    end if
    
    if (.not. s%piezometer .and. s%zOrd < 1) then
       write(*,*) 'ERROR: # of quadrature points at ',&
            & 'monitoring location must be > 0', s%zOrd
       stop
    end if

    if (s%timeseries) then
       ! if computing a time series, read time-related
       ! parameters from input file

       allocate(s%r(1),s%rD(1)) ! only one observation distance
       if (rval > w%rw) then
          s%r(1) = rval
       else
          write(*,*) 'ERROR: r must be > rw',rval
          stop
       end if
       
       if (s%piezometer) then
          s%zOrd = 1
       end if

       allocate(s%z(s%zOrd), s%zD(s%zOrd))       
       if (s%piezometer) then
          s%z(1) = s%zTop*f%b
       else
          if (s%zOrd > 1) then
             ! calc points spread out evenly along obs well screen
             s%z(1:s%zOrd) = linspace(s%zBot, s%zTop, s%zOrd)*f%b
          else
             ! one point goes to middle of interval
             s%z(1) = (s%zBot + s%zTop)*f%b/2.0 
          end if
       end if

       open(unit=22, file=trim(timeFileName), status='old', &
            & action='read',iostat=ioerr)
       if(ioerr /= 0) then
          write(*,'(A)') 'ERROR opening time input file '// &
               & trim(timeFileName)//' for reading'
       end if
       
       ! compute times?, # times to read (not used if vector computed)
       read(22,*) computeTimes, numTFile

       ! log_10(t_min), log_10(t_max), # times
       read(22,*) minLogT, maxLogT, numTComp   
       
       if (computeTimes) then
          s%nt = numTComp
       else
          s%nt = numTFile
       end if

       allocate(s%t(s%nt), s%tD(s%nt), h%sv(s%nt))

       if (computeTimes) then
          ! computing times 
          s%t = logspace(minLogT,maxLogT,numTComp)
       else
          ! times listed one per line
          do i=1,numtfile
             read(22,*,iostat=ioerr) s%t(i)
             if(ioerr /= 0) then
                write(*,*) 'ERROR reading time ',i,&
                     & ' from input file '//trim(timeFileName)
                stop
             end if
          end do
          close(22)
          if(any(s%t < spacing(0.0))) then
             write(*,*) 'ERROR t must be > 0:',s%t
             stop                
          end if
       end if
    else
       ! if computing a contour map or profile
       ! read space related parameters from other file

       allocate(s%t(1), s%tD(1), h%sv(1))
       if (tval > spacing(0.0)) then
          s%t(1) = tval
       else
          write(*,*) 'ERROR: t must be >0',tval
       end if

       ! piezometer doesn't mean anything
       ! when computing contours / profile 

       open(unit=23, file=trim(spaceFileName), status='old', &
            & action='read',iostat=ioerr)
       if(ioerr /= 0) then
          write(*,'(A)') 'ERROR opening space input file '// &
               & trim(spaceFileName)//' for reading'
       end if

       read(23,*) computeSpace, numRFile, numZFile
       read(23,*) minR, maxR, numRComp
       read(23,*) minZ, maxZ, numZComp

       if (computeSpace) then
          s%nr = numRComp
          s%nz = numZComp
       else
          s%nr = numRFile
          s%nz = numZFile
       end if
       
       allocate(s%z(s%nz), s%zD(s%nz), &
            &   s%r(s%nr), s%rD(s%nr))
      
       if (computeSpace) then
          s%r = linspace(minR,maxR,numRComp)
          s%z = linspace(minZ,maxZ,numZcomp)
       else
          ! radial distances listed on one line
          read(23,*,iostat=ioerr) s%r(1:s%nr)
          if(ioerr /= 0) then
             write(*,*) 'ERROR reading radial distance ',&
                  & 'from input file '//trim(spaceFileName)
             stop
          elseif(any(s%r < w%rw)) then
             write(*,*) 'ERROR r must be >= rw',s%r
             stop
          end if
          ! vertical depths listed on one line
          read(23,*,iostat=ioerr) s%z(1:s%nz)
          if(ioerr /= 0) then
             write(*,*) 'ERROR reading vertical depths ',&
                  & 'from input file '//trim(spaceFileName)
             stop
          elseif(any(s%z < 0.0) .or. any(s%z > f%b)) then
             write(*,*) 'ERROR z must be in range 0<=>b',s%z
             stop
          end if
       end if
    end if
    
    s%nt = size(s%t,1)
    s%nz = size(s%z,1)
    s%nr = size(s%r,1)

    ! determine which layer z point(s) are in
    allocate(s%zLay(size(s%z)))
    where(s%z(:) <= w%l)
       ! above well screen
       s%zLay = 1
    elsewhere
       where(s%z(:) < w%d)
          ! beside well screen
          s%zLay = 2
       elsewhere
          ! below well screen
          s%zLay = 3
       end where
    end where
    
    ! output filename
    read(19,*) s%outfilename                        
    close(19)

    ! ## compute dimensionless quantities #####

    ! characteristic length / time / head
    s%Lc = f%b
    s%Tc = f%b**2/(f%Kr/f%Ss)
    s%Hc = 2*PI*f%Kr*f%b/w%Q  
    
    ! compute derived or dimensionless properties
    f%sigma = f%Sy/(f%Ss*f%b)
    f%alphaD = f%kappa/f%sigma

    ! pumping and monitoring well screens are entered relative to
    ! thickness of aquifer (which may not be the characteristic length - see above)

    ! dimensionless lengths
    w%l = w%l*f%b
    w%d = w%d*f%b
    w%lD = w%l/s%Lc
    w%dD = w%d/s%Lc
    w%bD = abs(w%lD - w%dD)

    s%zD(:) = s%z(:)/s%Lc
    s%rD(:) = s%r(:)/s%Lc

    ! dimensionless times
    s%tD(:) = s%t(:)/s%Tc

    ! ## echo computed quantities #####

    if (.not. s%quiet) then
       write(*,'(A,'//RFMT//')') 'kappa:   ',f%kappa
       write(*,'(A,'//RFMT//')') 'sigma:   ',f%sigma
       write(*,'(A,'//RFMT//')') 'alpha_D: ',f%alphaD
       write(*,'(A,'//RFMT//')') 'Tc:  ',s%Tc
       write(*,'(A,'//RFMT//')') 'Lc:  ',s%Lc
       write(*,'(A,'//RFMT//')') 'b_D: ',w%bD
       write(*,'(A,'//RFMT//')') 'l_D: ',w%lD
       fmt = '(A,I0,1X,     (ES09.03,1X))           '
       write(fmt(10:14),'(I5.5)') s%nz
       write(*,fmt) 'z  : ',s%nz,s%z
       write(*,fmt) 'z_D: ',s%nz,s%zD
       write(fmt(10:14),'(I5.5)') s%nr
       write(*,fmt) 'r  : ',s%nr,s%r 
       write(*,fmt) 'r_D: ',s%nr,s%rD
       write(fmt(10:14),'(I5.5)') s%nt
       write(*,fmt) 't  : ',s%nt,s%t
       write(*,fmt) 't_D: ',s%nt,s%tD
       write(*,'(A,2('//RFMT//',1X),I0)') 'relative monitoring screen '//&
            & 'zTop/b, zBot/b ,zOrd: ', s%zTop, s%zBot, s%zOrd
       write(*,'(A,I0,2('//RFMT//',1X))') 'deHoog: M,alpha,tol: ', &
            & lap%M, lap%alpha, lap%tol
       write(*,'(A,2(I0,1X))'), 'tanh-sinh: k, num extrapolation steps ', &
            & ts%k, ts%nst
       write(*,'(A,4(I0,1X))'), 'GL: J0 split, num zeros accel, GL-order ',&
            & h%j0s(:), gl%nacc, gl%ord
    end if

    terms = maxval(h%j0s(:)) + gl%nacc + 1
    allocate(h%j0z(terms))

    ! ## compute zeros of J0 bessel function #####

    ! asymptotic estimate of zeros - initial guess
    forall (i=0:terms-1) h%j0z(i+1) = (i + 0.75)*PIEP
    do i=1,terms
       x = h%j0z(i)
       NR: do
          ! Newton-Raphson (f2008 bessel function name convention)
          dx = bessel_j0(x)/bessel_j1(x)
          x = x + dx
          if(abs(dx) < spacing(x)) then
             exit NR
          end if
       end do NR
       h%j0z(i) = x
    end do
    
    ! split between finite/infinite part should be 
    ! small for large time, large for small time
    zRange = maxval(h%j0s(:)) - minval(h%j0s(:))
    minLSp = floor(minval(log10(s%tD)))   ! min log(td) -> maximum split
    maxLSp = ceiling(maxval(log10(s%tD))) ! max log(td) -> minimum split
    spRange = maxlsp - minlsp + 1 
    h%sv = minval(h%j0s(:)) + int(zrange*((maxlsp - log10(s%tD))/spRange))

  end subroutine read_input

  subroutine write_timeseries_header(w,f,s,lap,h,gl,ts,unit)
    use constants, only : RFMT
    use types, only : invLaplace, invHankel, GaussLobatto, tanhSinh, well, formation, solution
    
    type(invLaplace), intent(in) :: lap
    type(invHankel), intent(in) :: h
    type(GaussLobatto), intent(in) :: gl
    type(TanhSinh), intent(in) :: ts
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    type(solution), intent(in) :: s
    integer, intent(in) :: unit
    integer :: ioerr
    character(32) :: fmt

    open (unit=unit, file=s%outFileName, status='replace', action='write', iostat=ioerr)
    if (ioerr /= 0) then
       write(*,'(A)') 'cannot open output file '//trim(s%outFileName)//' for writing' 
       stop
    end if
  
    ! echo input parameters at head of output file
    write(unit,'(A)') '# -*-auto-revert-*-'
    write(unit,'(A,3(L1,1X))') '# dimensionless?, timeseries?, piezometer? :: ', &
         & s%dimless, s%timeseries, s%piezometer
    write(unit,'(A,'//RFMT//')') '# Q (volumetric pumping rate) :: ', &
         & w%Q
    write(unit,'(A,'//RFMT//')') '# b (initial sat thickness) :: ', &
         & f%b
    write(unit,'(A,2('//RFMT//',1X))') '# l/b,d/b (screen bot & top) :: ',&
         & w%l, w%d  
    write(unit,'(A,2('//RFMT//',1X))') '# rw,rc (well/casing radii) :: ',&
         & w%rw, w%rc  
    write(unit,'(A,2('//RFMT//',1X))') '# Kr,kappa :: ', f%Kr, f%kappa
    write(unit,'(A,2('//RFMT//',1X))') '# Ss,Sy :: ',f%Ss, f%Sy
    write(unit,'(A,'//RFMT//')') '# gamma :: ',f%gamma
    fmt = '(A,I0,A,    ('//RFMT//',1X))       '
    write(fmt(9:12),'(I4.4)') size(lap%timePar)
    write(unit,fmt) '# pumping well time behavior :: ',lap%timeType, &
         & trim(lap%timeDescrip(lap%timeType)), lap%timePar
    write(unit,'(A,I0,2('//RFMT//',1X))') '# deHoog M, alpha, tol :: ',&
         & lap%M, lap%alpha, lap%tol
    write(unit,'(A,2(I0,1X))'), '# tanh-sinh: k, n extrapolation steps :: ',&
         & ts%k, ts%nst
    write(unit,'(A,4(I0,1X))'), '# GLquad: J0 split, n 0-accel, GL-order :: ',&
         & h%j0s(:), gl%nacc, gl%ord
    if(s%piezometer) then
       write(unit,'(A,2('//RFMT//',1X))') '# point obs well r,z :: ',s%r(1),s%z(1)
    else
       write(unit,'(A,3('//RFMT//',1X),I0)') '# screened obs well r,zTop,zBot,zOrd :: ',&
            & s%r(1), s%zTop, s%zBot, s%zOrd
    end if
    write(unit,'(A,I0)') '# times :: ',s%nt
    write(unit,'(A,2('//RFMT//',1X))') '# characteristic length, time :: ',s%Lc,s%Tc
    if (s%dimless) then
       write (unit,'(A,/,A,/,A)') '#','#       t_D              '//&
            & trim(s%modelDescrip(s%model)), &
            & '#----------------------------------------'
       ! TODO add header for derivative wrt logt
    else
       write(unit,'(A,1'//RFMT//')') '# characteristic head ::',s%Hc
       write (unit,'(A,/,A,/,A)') '#','#       t                '//&
            & trim(s%modelDescrip(s%model)), &
            & '#----------------------------------------'
       ! TODO add header for derivative wrt logt
    end if

  end subroutine write_timeseries_header

  subroutine write_contour_header(w,f,s,lap,h,gl,ts,unit)
    use constants, only : RFMT
    use types, only : invLaplace, invHankel, GaussLobatto, tanhSinh, well, formation, solution
    
    type(invLaplace), intent(in) :: lap
    type(invHankel), intent(in) :: h
    type(GaussLobatto), intent(in) :: gl
    type(TanhSinh), intent(in) :: ts
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    type(solution), intent(in) :: s
    integer, intent(in) :: unit
    integer :: ioerr
    character(32) :: fmt

    open (unit=unit, file=s%outFileName, status='replace', action='write', iostat=ioerr)
    if (ioerr /= 0) then
       write(*,'(A)') 'cannot open output file '//trim(s%outFileName)//' for writing' 
       stop
    end if
  
    ! echo input parameters at head of output file
    write(unit,'(A)') '# -*-auto-revert-*-'
    write(unit,'(A,2(L1,1X))') '# dimensionless?, timeseries? :: ', &
         & s%dimless, s%timeseries
    write(unit,'(A,'//RFMT//')') '# Q (volumetric pumping rate) :: ', &
         & w%Q
    write(unit,'(A,'//RFMT//')') '# b (initial sat thickness) :: ', &
         & f%b
    write(unit,'(A,2('//RFMT//',1X))') '# l,d (screen bot & top) :: ',&
         & w%l, w%d  
    write(unit,'(A,2('//RFMT//',1X))') '# rw,rc (well/casing radii) :: ',&
         & w%rw, w%rc  
    write(unit,'(A,2('//RFMT//',1X))') '# Kr,kappa :: ', f%Kr, f%kappa
    write(unit,'(A,2('//RFMT//',1X))') '# Ss,Sy :: ',f%Ss, f%Sy
    write(unit,'(A,'//RFMT//')') '# gamma :: ',f%gamma
    fmt = '(A,I0,A,    ('//RFMT//',1X))       '
    write(fmt(9:12),'(I4.4)') size(lap%timePar)
    write(unit,fmt) '# pumping well time behavior :: ',lap%timeType, &
         & trim(lap%timeDescrip(lap%timeType)), lap%timePar
    write(unit,'(A,I0,2('//RFMT//',1X))') '# deHoog M, alpha, tol :: ',&
         & lap%M, lap%alpha, lap%tol
    write(unit,'(A,2(I0,1X))'), '# tanh-sinh: k, n extrapolation steps :: ',&
         & ts%k, ts%nst
    write(unit,'(A,4(I0,1X))'), '# GLquad: J0 split, n 0-accel, GL-order :: ',&
         & h%j0s(:), gl%nacc, gl%ord
    fmt = '(A,I0,1X,    ('//RFMT//',1X))     '
    write(fmt(10:13),'(I4.4)') s%nr
    write(unit,fmt) '# num r locations, rlocs :: ',s%nr, s%r(:)
    write(fmt(10:13),'(I4.4)') s%nz
    write(unit,fmt) '# num z locations, zlocs :: ',s%nz, s%z(:)
    write(unit,'(A,'//RFMT//')') '# time :: ',s%t(1)
    if (s%dimless) then
       write (unit,'(A,/,A,/,A)') '#','#         z_D           r_D           '&
            & //'             '// trim(s%modelDescrip(s%model)), &
            & '#------------------------------------------------------------'
       ! TODO add header for derivative wrt logt
    else
       write (unit,'(A,/,A,/,A)') '#','#         z              r            '&
            & //'             '// trim(s%modelDescrip(s%model)), &
            & '#------------------------------------------------------------'
       ! TODO add header for derivative wrt logt
    end if

  end subroutine write_contour_header
end module driver_io

