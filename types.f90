
module types
  use constants, only : DP, EP
  implicit none

  public

  ! Inverse Laplace Transform parameters
  type :: invLaplace
     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha = -999., tol = -999.

     ! number of Fourier series terms
     integer :: M = -999

     ! length of solution vector (2*M+1)
     integer :: np = -999

     complex(EP), allocatable :: p(:)

     ! time behavior / parameters 
     ! 1 = step on,              tpar(1)  = on time
     ! 2 = finite pulse,         tpar(1:2) = on/off time
     ! 3 = instan. pulse         tpar(1) = pulse time
     ! 4 = stairs,               tpar(1) = time step (increasing Q by integer multiples 
     !                                     @ integer multiples tpar(1)); tpar(2) =off time.
     ! 5 = + only square wave,   tpar(1) = 1/2 period of wave; tpar(2) = start time
     ! 6 = cosine(tpar(1)*t),    tpar(1) = frequency multiplier; tpar(2) = start time
     ! 7 = + only tri wave,      tpar(1) = 1/4 period of wave; tpar(2) = start time
     ! 8 = +/- square wave,      tpar(1) = 1/2 period of wave; tpar(2) = start time
     ! n<0 = arbitrary piecewise constant rate, comprised of n steps from tpar(1) to tfinal
     !                tpar(1:n) = starting times of each step
     !                tpar(n+1) = final time of last step
     !                tpar(n+2:2*n+1) = strength at each of n steps 
     ! (is multiplied by constant strength too -- you probably want to set that to unity)
     character(80), parameter, dimension(9) :: timeDescrip = &
          & ['step on; tpar(1) = on time',&
          &  'finite pulse; tpar(1:2) = on/off time',&
          &  'infinitessimal pulse; tpar(1) = pulse location',&
          &  'stairs; tpar(1) = time step (Q increase by integer multiples); tpar(2) = off time',&
          &  'rectified square wave; tpar(1) = 1/2 period of wave; tpar(2) = start time',&
          &  'cos(omega*t); tpar(1) = omega; tpar(2) = start time',&
          &  'rectified triangular wave; tpar(1) = 1/4 period of wave; tpar(2) = start time',&
          &  'rectified square wave; tpar(1) = 1/2 period of wave; tpar(2) = start time',&
          &  'piecewise constant rate (n steps); tpar(1:n)=ti; tpar(n+1)=tfinal; tpar(n+2:)=Q']
     
     ! type of time behavior  (see above)
     integer :: timeType = -999

     ! parameters related to different time behaviors (on, off, etc)
     real(EP), allocatable :: timePar(:)

  end type invLaplace

  ! inverse Hankel transform parameters
  type :: invHankel
     ! zeros of J0 Bessel function
     real(DP), allocatable :: j0z(:) ! locations of zeros of J0 bessel fnc
     integer, allocatable :: sv(:) ! split index  vector

     ! min/max j0 split between infinite/fininte integrals
     integer, dimension(2) :: j0s = [-999, -999] 
  end type invHankel

  ! parameters specific to GL quadrature
  type :: GaussLobatto
     integer :: nacc = -999, err -999

     ! abcissa and weights
     real(EP), allocatable :: x(:), w(:)

     ! order of integration
     integer :: ord = -999
  end type GaussLobatto

  ! parameters specific to tanh-sinh quadrature  
  type :: TanhSinh
     integer :: k = -999, N = -999, nst = -999
     real(EP), allocatable :: w(:), a(:), hh(:)
     integer, allocatable :: kk(:), NN(:), ii(:)
     
     ! error in polynomial extrapolation
     complex(EP) :: polerr = (-999., -999.)
  end type TanhSinh

  ! parameters related to well/completion
  type :: well
     real(DP) :: l = -999. ! aquifer top to screen/packer top dist.
     real(DP) :: d = -999. ! aquifer top to screen/packer bottom dist.
     real(DP) :: rw = -999., rc = -999. ! well / casing radii
     real(DP) :: Q = -999. ! volumetric pumping rate
 
     ! dimensionless parameters
     real(DP) :: lD = -999.  ! dimensionless l
     real(DP) :: dD = -999.  ! dimensionless d
     real(DP) :: bD = -999.  ! dimensionless screen length
  end type well

  ! parameters related to formation/aquifer
  type :: formation
     real(DP) :: b = -999.  ! aquifer thickness
     real(DP) :: Kr = -999. ! radial hydraulic conductivity
     real(DP) :: kappa  = -999.  ! Kz/Kr ratio
     real(DP) :: Ss = -999. ! specific storage
     real(DP) :: Sy = -999. ! specific yield
     real(DP) :: gamma = -999.  ! dimensionless skin (1=no skin)
     real(DP) :: usL = -999. ! thickness of unsaturated zone
     real(DP) :: usalpha = -999. ! unzaturated zone sorbtive number

     ! computed aquifer parameters
     real(DP) :: sigma = -999.  ! Sy/(Ss*b)
     real(DP) :: alphaD = -999. ! kappa/sigma
  end type formation

  ! parameters related to numerical solution
  type(time) :: solution
     real(DP) :: Lc = -999.  ! characteristic length
     real(DP) :: Tc = -999.  ! characteristic time

     ! which unconfined model to use?
     integer :: model = -999
     ! 0 = Theis (confined fully penetrating)
     ! 1 = Hantush (confined partially penetrating)
     ! 2 = Boulton 195?
     ! 3 = Neuman 1974 
     ! 4 = Moench 199?
     ! 5 = Mishra/Neuman 2011
     ! 6 = Malama 2011
     
     character(13), parameter, dimension(0:6) :: modelDescrip = [ &
          & 'Theis', 'Hantush', 'Boulton', 'Neuman 74',&
          & 'Moench', 'Mishra/Neuman', 'Malama']

     logical :: quiet = .false.  ! output debugging to stdout?
     logical :: dimless = .false.  ! output dimensionless solution?
     logical :: timeseries = .false. ! vector of times, one location?
     logical :: piezometer = .false. ! point observation location?

     ! number of times / locations to be computed or read from file
     integer :: nt = -999, nr = -999, nz = -999  

     ! vector of times / locations to compute solution at
     real(DP), allocatable :: t(:), tD(:)
     real(DP), allocatable :: r(:), rD(:) 
     real(DP), allocatable :: z(:), zD(:)

     ! top/bot of monitoring well screen 
     real(DP) :: zTop = -999., zBot = -999. 
     ! order of quadrature at monitoring well screen
     integer :: zOrd = -999 

     integer, parameter :: NUMCHAR = 128
     character(NUMCHAR) :: outfilename, infilename

     character(7) :: rfmt = 'ES14.07'

  end type solution

end module types

