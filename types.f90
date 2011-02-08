
module types
  use constants, only : DP
  implicit none

  private
  public :: invlt, well, fm, soln

  type :: invLaplace
     ! Inverse Laplace Transform parameters

     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha = -999., tol = -999.

     ! number of Fourier series terms
     integer :: M = -999

     ! length of solution vector (2*M+1)
     integer :: np

  end type invLaplace

  type :: invHankel
     ! inverse Hankel transform parameters

     ! zeros of J0 Bessel function
     real(DP), allocatable :: j0zero(:)
     integer :: splitrange = -999, zerorange = -999
     integer, allocatable :: splitv(:)

     integer, dimension(2) :: j0split = [-999, -999]

  end type invHankel

  type :: GaussLobatto
     ! parameters specific to GL quadrature

     integer :: naccel = -999, err -999

     ! abcissa and weights
     real(EP), allocatable :: GLx(:),GLw(:)

     ! order of integration
     integer :: order = -999

     ! error in polynomial extrapolation
     complex(EP) :: polerr

  end type GaussLobatto
  
  type :: TanhSinh
     ! parameters specific to tanh-sinh quadrature

     integer :: k = -999, N = -999, nst = -999
     complex(EP), allocatable :: tmp(:,:)
     real(EP), allocatable :: w(:), a(:), hh(:)
     integer, allocatable :: kk(:), NN(:), ii(:)
     
  end type TanhSinh

  type :: well
     ! parameters related to well/completion
     real(DP) :: l = -999. ! aquifer top to screen/packer top dist.
     real(DP) :: d = -999. ! aquifer top to screen/packer bottom dist.
     real(DP) :: rw = -999., rc = -999. ! well / casing radii

     ! dimensionless parameters
     real(DP) :: lD = -999.  ! dimensionless l
     real(DP) :: dD = -999.  ! dimensionless d
     real(DP) :: bD = -999.  ! dimensionless screen length
     real(DP) :: rDw = -999. ! dimensionless rw

  end type well

  type :: formation
     ! parameters related to formation/aquifer
          
     real(DP) :: b = -999.  ! aquifer thickness
     real(DP) :: Kr = -999. ! radial hydraulic conductivity
     real(DP) :: kappa  = -999.  ! Kz/Kr ratio
     real(DP) :: Ss = -999. ! specific storage
     real(DP) :: Sy = -999. ! specific yield
     real(DP) :: gamma = -999.  ! dimensionless skin (1=no skin)

     ! computed aquifer parameters
     real(DP) :: sigma = -999.  ! Sy/(Ss*b)
     real(DP) :: alphaD = -999. ! kappa/sigma

  end type formation

  type :: solution
     ! parameters related to numerical solution

     real(DP) :: Lc = -999.  ! characteristic length
     real(DP) :: Tc = -999.  ! characteristic time

     logical :: quiet = .false.  ! output debugging to stdout?
     logical :: dimless = .false.  ! output dimensionless solution?

     ! either the number of times to be computed,
     ! or the number of times read from file
     integer :: nt = -999  

     ! vector of times to compute solution at
     real(DP), allocatable :: t(:), tD(:)

     integer, parameter :: NUMCHAR = 128
     character(NUMCHAR) :: outfilename, infilename, timefilename

     character(6) :: rfmt = 'ES14.7'

  end type solution

end module types

