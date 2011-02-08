module io
implicit none
private
public :: read_input, write_header

contains
  subroutine read_input(w,f,s,l,h,gl,ts)
    use constants, only : DP, PI
    use types
    use utilities, only : logspace
    
    type(invLaplace), intent(inout) :: l
    type(invHankel), intent(inout) :: h
    type(GaussLobatto), intent(inout) :: gl
    type(TanhSinh), intent(inout) :: ts
    type(well), intent(inout) :: w
    type(formation), intent(inout) :: f
    type(solution), intent(inout) :: s

    integer :: ioerr, terms
    integer :: zerorange, minlogsplit, maxlogsplit, splitrange
    integer :: numtfile, numtcomp, minlogt, maxlogt
    real(DP) :: x, dx

    intrinsic :: get_command_argument

    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! read input parameters from file some minimal error checking
    call get_command_argument(1,infile)
    if (len_trim(infile) == 0) then
       ! default input file name if no command-line argument
       infile = 'input.dat'
    end if

    open(unit=19, file=infile, status='old', action='read',iostat=ioerr)
    if(ioerr /= 0) then
       write(*,'(A)') 'ERROR opening main input file '//trim(infilename)//' for reading'
       stop
    end if
  
    ! suppress output to screen, write dimensionless output?
    read(19,*) s%quiet, s%dimless  

    ! initial saturated aquifer thickness [L]
    read(19,*) f%b
  
    ! distance from aquifer top to bottom & top of packer (l > d) [L]
    read(19,*) w%l, w%d

    ! well and casing radii [L]
    read(19,*) w%rw, w%rc         

    ! aquifer radial K [L^2/T] and anisotropy ratio (Kz/Kr)
    read(19,*) f%Kr, f%kappa
    
    ! aquifer specific storage [1/L] and specific yield [-]
    read(19,*) f%Ss, f%Sy

    ! dimensionless skin
    read(19,*) f%gamma

    if (.not. s%quiet) then
       write(*,'(A,2(L1,1X))') 'quiet?, dimless output? ', s%quiet, s%dimless
       write(*,'(A,'//s%rfmt//')') 'b (initial aquier sat thickness):',f%b  
       write(*,'(A,2('//s%rfmt//',1X))') 'l,d (screen bot&top from aquifer top):',w%l,w%d  
       write(*,'(A,2('//s%rfmt//',1X))') 'rw,rc (well and casing radii):',w%rw,w%rc  
       write(*,'(A,3('//s%rfmt//',1X))') 'Kr,kappa,gamma',f%Kr, f%kappa, f%gamma
       write(*,'(A,2('//s%rfmt//',1X))') 'Ss,Sy',f%Ss, f%Sy
    end if

    if(any([w%b,f%Kr,f%kappa,f%Ss,w%l,w%rw,w%rc] < spacing(0.0))) then
       write(*,'(A)') 'input error: zero or negative distance or aquifer parameters' 
       stop
    end if

    if (w%d >= w%l) then
       write(*,'(2(A,'//s%rfmt//'))') 'input error: l must be > d; l=',w%l,' d=',w%d
       stop
    end if
  
    ! Laplace transform (deHoog et al) parameters
    read(19,*) l%M, l%alpha, l%tol

    if (l%M < 2) then
       write(*,'(A,I0)')  'deHoog # FS terms must be >= 1 M=',l%M
       stop 
    end if

    if (l%tol < epsilon(l%tol)) then
       l%tol = epsilon(l%tol)
       write(*,'(A,'//s%rfmt//')') 'WARNING: increased INVLAP tolerance to ',l%tol 
    end if

    ! tanh-sinh quadrature parameters
    read(19,*) ts%k, ts%nst
    
    ! Gauss-Lobatto quadrature parameters
    read(19,*) h%j0split(1:2), gl%naccel, gl%order  

    if(ts%k - ts%nst < 2) then
       write(*,'(2(A,I0),A)') 'Tanh-Sinh k is too low (',ts%k,') for given level&
            & of Richardson extrapolation (',ts%nst,').  Increase tsk or decrease nst.'
       stop
    end if

    if(any([h%j0split(:),gl%naccel, ts%k] < 1)) then
       write(*,'(A,4(I0,1X))') 'max/min split, # accelerated terms, ' // &
            & 'and tanh-sinh k must be >= 1:',h%j0split(:),gl%naccel, ts%k
       stop
    end if

    terms = maxval(h%j0split(:)) + gl%naccel + 1
    allocate(h%j0zero(terms))

    ! log_10(t_min), log_10(t_max), # times
    read(19,*) minlogt, maxlogt, numtcomp   

    ! compute times?, # times in file, time file          
    read(19,*) s%computetimes, numtfile, s%timefilename    

    ! output filename
    read(19,*) s%outfilename                        
    
    if (s%computetimes) then
       s%numt = numtcomp
    else
       s%numt = numtfile
    end if

    ! solution vectors
!    allocate(finint(2*lap%M+1), infint(2*lap%M+1), totlap(2*lap%M+1), p(2*lap%M+1),&
!         & splitv(numt), dt(numt))

    allocate(s%t(s%numt), s%tD(s%numt), splitv(s%numt))

    if (s%computetimes) then
       ! computing times 
       s%t = logspace(minlogt,maxlogt,numtcomp)
    else
       ! read times from file
       open(unit=22, file=s%timefilename, status='old', action='read',iostat=ioerr)
       if(ioerr /= 0) then
          write(*,'(A)') 'ERROR opening time data input file '//trim(timefilename) &
               & //' for reading'
          stop 
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
       end if
       close(22)
    end if

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
       write(*,'(A,I0,2('//s%rfmt//',1X))') 'deHoog: M,alpha,tol: ',l%M, l%alpha, l%tol
       write(*,'(A,2(I0,1X))'), 'tanh-sinh: k, num extrapolation steps ',ts%k, ts%nst
       write(*,'(A,4(I0,1X))'), 'GL: J0 split, num zeros accel, GL-order ',h%j0split(:),&
            & gl%naccel, gl%order
       write(*,'(A,L1)') 'compute times? ',s%computetimes
       if(s%computetimes) then
          write(*,'(A,3(I0,1X))') 'log10(tmin), log10(tmax), num times ',&
               & minlogt,maxlogt,numtcomp
       else
          write(*,'(A,I0,1X,A)') 'num times, filename for t inputs ',&
               & numtfile,trim(timefname)
       end if
    end if

    ! compute zeros of J0 bessel function
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
 
    !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! finite portion of Hankel integral for each time level
    
    l%np = 2*l%M+1  ! number of Laplace transform Fourier series coefficients

    ts%N = 2**ts%k - 1
    allocate(ts%w(ts%N),ts%a(ts%N),fa(ts%N,l%np),ii(ts%N),tmp(ts%nst,l%np),&
         & ts%kk(ts%nst),ts%NN(ts%nst),ts%hh(ts%nst))

  end subroutine read_input

  subroutine write_header(w,f,s,l,h,gl,ts,unit)
    use types
    
    type(invLaplace), intent(in) :: l
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
    write(20,'(A,L1)') '# quiet? ', s%quiet
    write(20,'(A,'//s%rfmt//')') '# b (initial aquier sat thickness): ',f%b
    write(20,'(A,2('//s%rfmt//',1X))') '# l,d (packer bot & top from top of aquifer): ',w%l,w%d  
    write(20,'(A,2('//s%rfmt//',1X))') '# rw,rc (well and casing radii): ',w%rw,w%rc  
    write(20,'(A,2('//s%rfmt//',1X))') '# Kr,kappa: ',f%Kr, f%kappa
    write(20,'(A,2('//s%rfmt//',1X))') '# Ss,Sy: ',f%Ss, f%Sy
    write(20,'(A,I0,2('//s%rfmt//',1X))') '# deHoog M, alpha, tol: ',l%M, l%alpha, l%tol
    write(20,'(A,2(I0,1X))'), '# tanh-sinh: k, num extrapolation steps: ',ts%k, ts%nst
    write(20,'(A,4(I0,1X))'), '# GL quad: J0 split, # zeros accel, GL-order: ',h%j0split(:), gl%naccel, gl%order
    write(20,'(A,L1)') '# compute times?: ',s%computetimes
    if (s%dimless) then
       write (20,'(A,/,A,/,A)') '#','#         t_D                       H^{-1}[ L^{-1}[ f_D ]] ',&
            & '#------------------------------------------------------------'
    else
       write (20,'(A,/,A,/,A)') '#','#         t                         H^{-1}[ L^{-1}[ f ]]   ',&
            & '#------------------------------------------------------------'
    end if

  end subroutine write_header
end module io

