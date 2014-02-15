!
! Copyright (c) 2012 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module driver_io
implicit none

private
public :: read_input, write_timeseries_header, write_contour_header

contains
  subroutine read_input(w,f,s,lap,h,gl,ts)
    use constants, only : EP, DP, PI, PIEP, NUMCHAR, RFMT, SFMT
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
    character(43) :: fmt

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

    open(unit=19, file=inputFileName, status='old', action='read', iostat=ioerr)
    if(ioerr /= 0) then
       write(*,*) 'ERROR opening main input file '//trim(inputFileName)//&
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

    if(s%model < 0 .or. s%model > 6) then
       write(*,'(A,I0,A)') 'ERROR invalid model choice ',s%model,' valid models are:'
       do i=0,6
          write(*,'(I0,1X,A)') i,trim(s%modelDescrip(i))
       end do
       write(*,'(A)') 'Malama models with beta=0 correspond to Neuman 72/74'
       stop
    end if

    ! ## pumping well parameters #####

    ! volumetric pumping rate [L^3/T]
    read(19,*) w%Q

    ! __distance_from_aquifer_top__ to bottom & top of packer (l > d) [L]
    ! NOT THE Z-COORDINATE OF THE TOP OR BOTTOM OF SCREEN!
    read(19,*) w%l, w%d

    ! pumping well: well and casing radii [L]
    read(19,*) w%rw, w%rc

    ! dimensionless skin
    read(19,*) f%gamma

    ! pumping well time behavior
    read(19,*) lap%timeType
    backspace(19) ! re-read this record after allocating space
    if (lap%timeType > -1) then
       ! functional time behavior (two or fewer parameters)
       allocate(lap%timePar(2))
       lap%timePar = -999.
       read(19,*) lap%timeType,lap%timePar(:)
    else
       ! arbitrarily long piecewise-constant time behavior
       allocate(lap%timePar(-2*lap%timeType+1))
       lap%timePar = -999.
       read(19,*) lap%timeType,lap%timePar(:)
    end if

    ! ## aquifer / formation parameters #####

    ! initial saturated aquifer thickness [L]
    read(19,*) f%b

    ! aquifer radial K [L^2/T] and anisotropy ratio (Kz/Kr)
    read(19,*) f%Kr, f%kappa

    ! aquifer specific storage [1/L] and specific yield [-]
    read(19,*) f%Ss, f%Sy

    ! Malama linearization beta parameter
    ! Moench M (number of alpha coefficients)
    read(19,*) f%beta, f%MoenchAlphaM
    backspace(19)
    if (f%MoenchAlphaM < 0) then
       write(*,*) 'ERROR: number of Moench alpha parameters must be >1',f%MoenchAlphaM
    end if
    allocate(f%MoenchAlpha(f%MoenchAlphaM),f%MoenchGamma(f%MoenchAlphaM))
    f%MoenchAlpha = -999.
    f%MoenchGamma = -999.
    read(19,*) f%beta, f%MoenchAlphaM, f%MoenchAlpha(:)

    ! Mishra/Neuman unsaturated model parameters
    ! capacity & conductivity sorptive numbers (1/length)
    ! negative air-entry & saturation pressures (>0)
    ! unsaturated zone thickness (length)
    read(19,*) f%ac, f%ak, f%psia, f%psik, f%usL, s%MNtype, s%order

    ! ## account for differences between Malama's implementation and original ##
    ! reset variables here (even though not used), so the writing of these to 
    ! output file headers and echoing input to screen properly shows
    ! what they are effectively assumed to be.
    ! (so you hopefully can notice this even with s%quiet = 0)
    if (s%MNtype == 1) then
       if (abs(f%ac - f%ak) > epsilon(1.0)) then
          if (s%quiet > 0) then
             write(*,'(A)') 'WARNING1: Mishra-Neuman implementation '//&
                  & '1 assumes ac=ak: using ak (ignoring ac)'
             write(*,'(2(A,'//RFMT//'))') 'WARNING1: you gave ac=',f%ac,' ak=',f%ak
          endif
          f%ac = f%ak
       end if
       if (abs(w%l - f%b) > epsilon(1.0)) then
          if (s%quiet > 0) then
             write(*,'(A)') 'WARNING2: Mishra-Neuman implementation '//&
                  & '1 assumes pumping well screen goes to bottom of formation (l=b)'
             write(*,'(2(A,'//RFMT//'))') 'WARNING2: you gave l=',w%l,' b=',f%b
          endif
          w%l = f%b
       end if
       if (w%d > epsilon(1.0)) then
          if (s%quiet > 0) then
             write(*,'(A)') 'WARNING3: Mishra-Neuman implementation '//&
                  & '1 assumes pumping well screen goes to top of formation (d=0)'
             write(*,'(A,'//RFMT//')') 'WARNING3: you gave d=',w%d
          endif
          w%d = 0.0
       end if
    end if
    
    ! ## echo check parameters #####

    if (s%quiet > 1) then
       write(*,'(A,A,1X,L1)') 'model, dimless output?:: ',&
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
          write(fmt(8:11), '(I4.4)') size(lap%timePar(:-lap%timeType+1),1)
          write(fmt(26:29),'(I4.4)') size(lap%timePar( -lap%timeType+2:),1)
          write(*,'(A)') 'pumping well time behavior :: '//trim(lap%timeDescrip(9))
          write(*,fmt) 'time behavior:  ti, tf | Q multiplier each step :: ', &
               & lap%timeType,lap%timePar(:-lap%timeType+1),'| ',&
               & lap%timePar(-lap%timeType+2:)
       end if
       write(*,'(A,'//RFMT//')') 'b (initial aquier sat thickness):: ',f%b
       write(*,'(A,2('//RFMT//',1X))') 'l,d (screen bot&top measured from above):: ', w%l, w%d
       write(*,'(A,2('//RFMT//',1X))') 'rw,rc (well and casing radii):: ', w%rw, w%rc
       write(*,'(A,3('//RFMT//',1X))') 'Kr,kappa,gamma:: ', f%Kr, f%kappa, f%gamma
       write(*,'(A,2('//RFMT//',1X))') 'Ss,Sy:: ', f%Ss, f%Sy
       if (s%model == 3) then
          fmt = '(A,I0,    (1X,'//RFMT//'))                  '
          write(fmt(7:10),'(I4.4)') f%MoenchAlphaM
          write(*,fmt) 'alpha (Moench Delayed Yield):: ',f%MoenchAlphaM,f%MoenchAlpha(:)
       elseif (s%model == 4 .or. s%model == 5) then
          write(*,'(A,'//RFMT//')') 'beta (Malama linearization factor):: ',f%beta
       elseif (s%model == 6) then
          write(*,'(A,4('//RFMT//',1X))') 'Mishra&Neuman ac,ak (sorptive #s), psia,psik '//&
               & '(neg air-entry & sat. pressures):: ', f%ac,f%ak,f%psia,f%psik
          write(*,'(A,I0)') 'Mishra&Neuman type of solution (0=naive,2=FD)::',s%MNtype
          write(*,'(A,'//RFMT//')') 'unsaturated zone thickness:: ',f%usL
          write(*,'(A,I0,'//RFMT//')') 'unsaturated zone FD order, FD h:: ',&
               & s%order,f%usL/(s%order-1)
       end if
    end if

    if(s%model > 0 .and. any([f%gamma, w%d, w%l] < 0.0)) then
       write(*,*) 'ERROR: negative geometry parameters:', &
            &  f%gamma, w%d, w%l
       stop
    end if

    if (any([f%b,f%Kr,f%Ss] < spacing(0.0))) then
       write(*,*) 'ERROR: zero or negative aquifer parametrs:',[f%b,f%Kr,f%Ss]
       stop
    end if

    if(s%model > 2 .and. any([f%kappa,f%Sy] < spacing(0.0))) then
       write(*,*) 'ERROR: zero or negative unconfined aquifer parameters:', &
            &  [f%kappa,f%Sy]
       stop
    end if

    if (s%model > 0 .and. w%d >= w%l) then
       write(*,*) 'ERROR: screen top/bottom l must be > d; l=',w%l,' d=',w%d
       stop
    end if

    if (s%model == 6) then
       if (any([f%ac,f%ak,f%usL,f%psia,f%psik] < 0.0)) then
          write(*,*) 'ERROR: ivalid Mishra/Neuman parameters',&
               & f%ac,f%ak,f%usL,f%psia,f%psik
          stop
       end if
       if (s%order < 3) then
          write(*,*) 'ERROR: order of Mishra/Neuman finite difference matrix must be >=3',s%order
          stop
       end if

       if (s%MNtype < 0 .and. s%MNtype > 2) then
          write(*,*) 'ERROR: invalid choice for Mishra/Neuman solution type '&
               & //'(naive=0,Malama=1,finite difference=2)',s%MNtype
          stop
       end if
    end if

    if ((s%model == 4 .or. s%model == 5) .and. f%beta < 0.0) then
       write(*,*) 'ERROR: Malama linearization beta cannot be negative',f%beta
       stop
    end if

    if (s%model == 3 .and. any(f%MoenchAlpha(:) < 0.0)) then
       write(*,*) 'ERROR: Moench alpha parameters cannot be negative',f%MoenchAlpha
       stop
    end if


    ! ## numerical implementation-related parameters #####

    ! Laplace transform (deHoog et al) parameters
    read(19,*) lap%M, lap%alpha, lap%tol

    ! tanh-sinh quadrature parameters
    ! integration order 2^(k-1) and Richardson extrapolation order
    read(19,*) ts%k, ts%R

    ! Gauss-Lobatto quadrature parameters
    ! max/min J0 zero to split at, # zeros to accelerate, GL-order
    read(19,*) h%j0s(1:2), gl%nacc, gl%ord

    ! ## checking of numerical parameters #####

    if (lap%M < 2) then
       write(*,*)  'ERROR: deHoog # FS terms must be >= 1 M=',lap%M
       stop
    end if

    if (lap%tol < epsilon(lap%tol)) then
       lap%tol = epsilon(lap%tol)
       write(*,'(A,'//RFMT//')') 'WARNING: increased INVLAP tolerance to ',&
            & lap%tol
    end if

    if(ts%k - ts%R < 2) then
       write(*,'(2(A,I0),A)') 'ERROR: Tanh-Sinh k is too low (',ts%k,&
            & ') for given level of Richardson extrapolation (',ts%R,&
            &').  Increase k or decrease nst.'
       stop
    end if
    if(ts%R < 1) then
       write(*,'(A,I0)') 'ERROR: Richardson extrapolation '//&
            &'level must be >= 1:  ',ts%R
       stop
    end if

    if(any([h%j0s(:),gl%nacc, ts%k] < 1)) then
       write(*,*) 'ERROR max/min split, # accelerated terms, ',&
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

    ! if computing contour or profile, these aren't used at all.
    ! if piezometer, midpoint between ztop and zbot is used.
    ! z is positive up, with zero at bottom of aquifer
    read(19,*) s%zTop, s%zBot, s%zOrd, s%rwobs, s%sF

    if (s%timeseries) then

       ! checking only related to timeseries (observation well parameters)
       if (s%zTop < s%zBot) then
          write(*,*) 'ERROR: for screened observation wells top of monitoring well'// &
               &'screen must be at or above bottom',s%zTop,s%zBot
          stop 666
       end if

       if (s%zTop > f%b .or. s%zBot < 0.0) then
          write(*,*) 'ERROR: top of monitoring well screen must be ',&
               & 'above bottom and both between 0 and b',s%zTop,s%zBot,f%b
          stop 667
       end if

       if (.not. s%piezometer .and. s%zOrd < 1) then
          write(*,*) 'ERROR: # of quadrature points at ',&
               & 'monitoring location must be > 0', s%zOrd
          stop
       end if

       if (s%rwobs < spacing(1.0)) then
          write(*,*) 'ERROR: monitoring well radius must be >0',s%rwobs
          stop 668
       end if
       if (s%sF < spacing(1.0)) then
          write(*,*) 'ERROR: monitoring well shape factor must be >0',s%sF
          stop 669
       end if

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
       s%zD = -999.
       ! calc points spread out evenly along obs well screen
       ! a single point is located in the middle of the wellscreen
       s%z(1:s%zOrd) = linspace(s%zBot, s%zTop, s%zOrd)

       open(unit=22, file=trim(timeFileName), status='old', &
            & action='read',iostat=ioerr)
       if(ioerr /= 0) then
          write(*,*) 'ERROR opening time input file '// &
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
       s%tD = -999.
       h%sv = -999

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
       s%tD = -999.
       h%sv = -999
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
          write(*,*) 'ERROR opening space input file '// &
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

       allocate(s%z(s%nz), s%zD(s%nz),s%r(s%nr), s%rD(s%nr))
       s%z = -999.
       s%zD = -999.
       s%r = -999.
       s%rD = -999.

       if (computeSpace) then
          s%r = linspace(minR,maxR,numRComp)
          if(minZ < 0.0 .or. maxZ > f%b) then
             write(*,*) 'ERROR z must be in range 0<=>b',minZ,maxZ
             stop
          end if
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

    ! output filename
    read(19,*) s%outfilename
    close(19)

    ! ## compute dimensionless quantities #####

    ! characteristic length / time / head
    s%Lc = f%b
    s%Tc = s%Lc**2/(f%Kr/f%Ss)
    s%Hc = w%Q/(4*PI*f%Kr*f%b)

    ! compute derived or dimensionless properties
    f%sigma = f%Sy/(f%Ss*f%b) ! this is the inverse of Neuman's (1972 & 1974) "sigma"
    f%alphaD = f%kappa/f%sigma
    f%betaD = f%beta/s%Lc

    ! dimensionless lengths
    w%lD = w%l/s%Lc
    w%dD = w%d/s%Lc
    w%bD = w%lD - w%dD
    w%rDw = w%rw/s%Lc
    s%rDwobs = s%rwobs/s%Lc

    f%MoenchGamma(:) = f%MoenchAlpha(:)*s%Tc

    ! dimensionless sorptive numbers
    f%acD = f%ac*s%Lc
    f%akD = f%ak*s%Lc
    f%lambdaD = (f%ak - f%ac)*s%Lc

    ! dimensionless pressures
    f%psiaD = f%psia/s%Lc
    f%psikD = f%psik/s%Lc
    f%usLD = f%usL/s%Lc
    f%b1 = f%psia - f%psik
    f%PsiD = f%b1/s%Lc

    ! dimensionlses coordinates
    s%zD(1:s%nz) = s%z(:)/s%Lc
    s%rD(1:s%nr) = s%r(:)/s%Lc
    s%tD(1:s%nt) = s%t(:)/s%Tc

    ! determine in which layer z point(s) are located
    ! z is positive up, with 0 at bottom
    ! l and d are positive down, with 0 at top
    allocate(s%zLay(size(s%zD)))
    s%zLay = -999

    where(s%zD(:) < spacing(1.0) .or. s%zD(:) < (1.0-w%lD))
       ! below well screen (or at bottom of aquifer zD==0)
       s%zLay = 1
    elsewhere
       where((s%zD(:)-1.0) > spacing(1.0) .or. s%zD(:) < (1.0-w%dD))
          ! beside/next-to well screen
          s%zLay = 2
       elsewhere
          ! above well screen (or top of aqufier)
          s%zLay = 3
       end where
    end where

    ! ## echo computed quantities #####

    if (s%quiet > 0) then
       write(*,'(A,'//RFMT//')') 'kappa   (Kz/Kr):   ',f%kappa
       write(*,'(A,'//RFMT//')') 'sigma   (Sy/S):   ',f%sigma
       write(*,'(A,'//RFMT//')') 'alpha_D:(kappa/sigma): ',f%alphaD
       write(*,'(A,2('//RFMT//',1X))') 'beta,beta_D: ',f%beta,f%betaD
       write(*,'(3(A,'//RFMT//'))') 'Tc:',s%Tc,' Lc:',s%Lc,' Hc:',s%Hc
       write(*,'(3(A,'//RFMT//'))') 'b_D:',w%bD,' l_D:',w%lD,' d_D:',w%dD
       fmt = '(A,I0,1X,     ('//SFMT//',1X))          '
       write(fmt(10:14),'(I5.5)') s%nz
       write(*,fmt) 'z  : ',s%nz,s%z
       write(*,fmt) 'z_D: ',s%nz,s%zD
       fmt = '(A,     (I0,1X))                       '
       write(fmt(4:8),'(I5.5)') s%nz
       write(*,fmt) 'zLay: ',s%zLay
       fmt = '(A,I0,1X,     ('//SFMT//',1X))          '
       write(fmt(10:14),'(I5.5)') s%nr
       write(*,fmt) 'r  : ',s%nr,s%r
       write(*,fmt) 'r_D: ',s%nr,s%rD
       write(fmt(10:14),'(I5.5)') s%nt
       write(*,fmt) 't  : ',s%nt,s%t
       write(*,fmt) 't_D: ',s%nt,s%tD
       write(*,'(A,2('//RFMT//',1X),I0)') 'monitoring screen '//&
            & 'zTop, zBot ,zOrd: ', s%zTop, s%zBot, s%zOrd
       write(*,'(A,2('//RFMT//',1X))') 'monitoring well rW and shape factor: ', &
            & s%rwobs, s%sF
       write(*,'(A,I0,2('//RFMT//',1X))') 'deHoog: M,alpha,tol: ', &
            & lap%M, lap%alpha, lap%tol
       write(*,'(A,2(I0,1X))') 'tanh-sinh: k, num extrapolation steps ', &
            & ts%k, ts%R
       write(*,'(A,4(I0,1X))') 'GL: J0 split, num zeros accel, GL-order ',&
            & h%j0s(:), gl%nacc, gl%ord
       write(*,'(A,4('//RFMT//',1X))') 'Mishra&Neuman acD,akD, psiaD,psikD :: ', &
            & f%acD,f%akD,f%psiaD,f%psikD
       write(*,'(A,3('//RFMT//',1X),2(I0,1X))') &
            & 'Mishra&Neuman LD, lambdaD, b1, MN type {0,1,2}, MN finite difference order:: ',&
            & f%usLD, f%lambdaD, f%b1, s%MNtype, s%order
    end if

    terms = maxval(h%j0s(:)) + gl%nacc + 1
    allocate(h%j0z(terms))

    ! ## compute zeros of J0 bessel function #####
    ! asymptotic estimate of zeros - initial guess
    forall (i=0:terms-1)
       h%j0z(i+1) = (i + 0.75)*PIEP
    end forall
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

  subroutine write_timeseries_header(w,f,s,lap,h,gl,ts,U)
    use constants, only : RFMT, EP
    use types, only : invLaplace, invHankel, GaussLobatto, tanhSinh, well, formation, solution

    type(invLaplace), intent(in) :: lap
    type(invHankel), intent(in) :: h
    type(GaussLobatto), intent(in) :: gl
    type(TanhSinh), intent(in) :: ts
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    type(solution), intent(in) :: s
    integer, intent(in) :: U
    integer :: ioerr
    character(32) :: fmt

    open (unit=U, file=s%outFileName, status='replace', action='write', iostat=ioerr)
    if (ioerr /= 0) then
       write(*,'(A)') 'cannot open output file '//trim(s%outFileName)//' for writing'
       stop
    end if

    ! echo input parameters at head of output file
    write(U,'(A)') '# -*-auto-revert-*-'
    write(U,'(A,I0,1X,A,I0)') '# model, EP :: ',s%model,trim(s%modelDescrip(s%model))//', ',EP
    write(U,'(A,3(L1,1X))') '# dimensionless?, timeseries?, piezometer? :: ', &
         & s%dimless, s%timeseries, s%piezometer
    write(U,'(A,'//RFMT//')') '# Q (volumetric pumping rate) :: ', &
         & w%Q
    write(U,'(A,'//RFMT//')') '# b (initial sat thickness) :: ', &
         & f%b
    write(U,'(A,2('//RFMT//',1X))') '# l,d (screen bot & top) :: ',&
         & w%l, w%d
    write(U,'(A,2('//RFMT//',1X))') '# rw,rc (well/casing radii) :: ',&
         & w%rw, w%rc
    write(U,'(A,2('//RFMT//',1X))') '# Kr,kappa :: ', f%Kr, f%kappa
    write(U,'(A,2('//RFMT//',1X))') '# Ss,Sy :: ',f%Ss, f%Sy
    write(U,'(A,'//RFMT//')') '# gamma :: ',f%gamma
    fmt = '(A,I0,A,    ('//RFMT//',1X))     '
    write(fmt(9:12),'(I4.4)') size(lap%timePar)
    write(U,fmt) '# pumping well time behavior :: ',lap%timeType, &
         & trim(lap%timeDescrip(lap%timeType)), lap%timePar
    write(U,'(A,I0,2('//RFMT//',1X))') '# deHoog M, alpha, tol :: ',&
         & lap%M, lap%alpha, lap%tol
    write(U,'(A,2(I0,1X))') '# tanh-sinh: k, n extrapolation steps :: ',&
         & ts%k, ts%R
    write(U,'(A,4(I0,1X))') '# GLquad: J0 split, n 0-accel, GL-order :: ',&
         & h%j0s(:), gl%nacc, gl%ord
    if(s%piezometer) then
       write(U,'(A,4('//RFMT//',1X))') '# point obs well r,rD,z,zD :: ',&
            &s%r(1),s%rD(1),s%z(1),s%zD(1)
    else
       write(U,'(A,3('//RFMT//',1X),I0)') '# screened obs well r,zTop,zBot,zOrd :: ',&
            & s%r(1), s%zTop, s%zBot, s%zOrd
       write(U,'(A,2('//RFMT//',1X))') '# screened obs well rW, shape factor :: ',&
            & s%rwobs, s%sF
    end if
    if(s%model == 4 .or. s%model == 5) then
       ! malama solutions
       write(U,'(A,'//RFMT//')') '# Malama beta linearization parameter :: ',f%beta
    elseif(s%model == 3) then
       fmt = '(A,I0,    (1X,'//RFMT//'))       '
       write(fmt(7:10),'(I4.4)') f%MoenchAlphaM
       write(*,fmt) '# Moench Delayed Yield (alpha):: ',f%MoenchAlphaM,f%MoenchAlpha(:)
    elseif(s%model == 6) then
       ! mishra/neuman solution
       write(U,'(A,5('//RFMT//',1X),I0)') '# Mishra/Neuman ac,ak,psia,psik,b1 ::',&
            & f%ac,f%ak,f%psia,f%psik,f%b1
       if (s%MNtype == 2) then
          write(U,'(A,I0,1X,'//RFMT//')') '# Mishra/Neuman FD order,h ::',&
               & s%order,f%usL/(s%order-1)
       elseif (s%MNtype == 1) then
          write(U, '(A'//RFMT//')') "# NB: Malama's Mishra/Neuman implementation (1) "//&
               &"assumes ac=ak and fully penetrating pumping well without wellbore storage"
       end if

    end if

    write(U,'(A,I0)') '# times :: ',s%nt
    write(U,'(A,2('//RFMT//',1X))') '# characteristic length, time :: ',s%Lc,s%Tc
    if (s%dimless) then
       write (U,'(A,/,A,/,A)') '#','#     t_D          '//&
            & trim(s%modelDescrip(s%model))//'         t*dh/d(log(t))', &
            & '#-------------------------------------------------------------'
    else
       write(U,'(A,1'//RFMT//')') '# characteristic head ::',s%Hc
       write (U,'(A,/,A,/,A)') '#','#     t            '//&
            & trim(s%modelDescrip(s%model))//'         t*dh/d(log(t))', &
            & '#-------------------------------------------------------------'
    end if

  end subroutine write_timeseries_header

  subroutine write_contour_header(w,f,s,lap,h,gl,ts,U)
    use constants, only : RFMT, EP
    use types, only : invLaplace, invHankel, GaussLobatto, tanhSinh, well, formation, solution

    type(invLaplace), intent(in) :: lap
    type(invHankel), intent(in) :: h
    type(GaussLobatto), intent(in) :: gl
    type(TanhSinh), intent(in) :: ts
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    type(solution), intent(in) :: s
    integer, intent(in) :: U
    integer :: ioerr
    character(32) :: fmt

    open (unit=U, file=s%outFileName, status='replace', action='write', iostat=ioerr)
    if (ioerr /= 0) then
       write(*,'(A)') 'cannot open output file '//trim(s%outFileName)//' for writing'
       stop
    end if

    ! echo input parameters at head of output file
    write(U,'(A)') '# -*-auto-revert-*-'
    write(U,'(A,I0,1X,A,I0)') '# model, EP :: ',s%model,trim(s%modelDescrip(s%model))//', ',EP
    write(U,'(A,2(L1,1X))') '# dimensionless?, timeseries? :: ', &
         & s%dimless, s%timeseries
    write(U,'(A,'//RFMT//')') '# Q (volumetric pumping rate) :: ', &
         & w%Q
    write(U,'(A,'//RFMT//')') '# b (initial sat thickness) :: ', &
         & f%b
    write(U,'(A,2('//RFMT//',1X))') '# l,d (screen bot & top) :: ',&
         & w%l, w%d
    write(U,'(A,2('//RFMT//',1X))') '# rw,rc (well/casing radii) :: ',&
         & w%rw, w%rc
    write(U,'(A,2('//RFMT//',1X))') '# Kr,kappa :: ', f%Kr, f%kappa
    write(U,'(A,2('//RFMT//',1X))') '# Ss,Sy :: ',f%Ss, f%Sy
    write(U,'(A,'//RFMT//')') '# gamma :: ',f%gamma
    fmt = '(A,I0,A,    ('//RFMT//',1X))     '
    write(fmt(9:12),'(I4.4)') size(lap%timePar)
    write(U,fmt) '# pumping well time behavior :: ',lap%timeType, &
         & trim(lap%timeDescrip(lap%timeType)), lap%timePar
    write(U,'(A,I0,2('//RFMT//',1X))') '# deHoog M, alpha, tol :: ',&
         & lap%M, lap%alpha, lap%tol
    write(U,'(A,2(I0,1X))') '# tanh-sinh: k, n extrapolation steps :: ',&
         & ts%k, ts%R
    write(U,'(A,4(I0,1X))') '# GLquad: J0 split, n 0-accel, GL-order :: ',&
         & h%j0s(:), gl%nacc, gl%ord
    fmt = '(A,I0,1X,    ('//RFMT//',1X))    '
    write(fmt(10:13),'(I4.4)') s%nr
    write(U,fmt) '# num r locations, rlocs :: ',s%nr, s%r(:)
    write(fmt(10:13),'(I4.4)') s%nz
    write(U,fmt) '# num z locations, zlocs :: ',s%nz, s%z(:)
    write(U,'(A,2('//RFMT//',1X))') '# time, tD :: ',s%t(1),s%tD(1)

    if(s%model == 4 .or. s%model == 5) then
       ! malama solutions
       write(U,'(A,'//RFMT//')') '# Malama beta linearization parameter :: ',f%beta
    elseif(s%model == 3) then
       fmt = '(A,I0,    (1X,'//RFMT//'))       '
       write(fmt(7:10),'(I4.4)') f%MoenchAlphaM
       write(*,fmt) '# Moench Delayed Yield (alpha):: ',f%MoenchAlphaM,f%MoenchAlpha(:)
    elseif(s%model == 6) then
       ! mishra/neuman solution
       write(U,'(A,5('//RFMT//',1X),I0)') '# Mishra/Neuman ac,ak,psia,psik,b1 ::',&
            & f%ac,f%ak,f%psia,f%psik,f%b1
       if (s%MNtype == 2) then
          write(U,'(A,I0,1X,'//RFMT//')') '# Mishra/Neuman FD order,h ::',&
               & s%order,f%usL/(s%order-1)
       elseif (s%MNtype == 1) then
          write(U, '(A'//RFMT//')') "# NB: Malama's Mishra/Neuman implementation (1) '//&
               &'assumes ac=ak and fully penetrating pumping well without wellbore storage"
       end if
    end if

    if (s%dimless) then
       write (U,'(A,/,A,/,A)') '#','#     z_D           r_D      '&
            & //'     '// trim(s%modelDescrip(s%model))//'          t*dh/d(log(t))', &
            & '#----------------------------------------------------------------------------'

    else
       write (U,'(A,/,A,/,A)') '#','#      z            r        '&
            & //'     '// trim(s%modelDescrip(s%model))//'          t*dh/d(log(t))', &
            & '#----------------------------------------------------------------------------'
    end if

  end subroutine write_contour_header
end module driver_io

