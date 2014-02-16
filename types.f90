!
! Copyright (c) 2012-2014 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module types
  use constants, only : DP, EP, NUMCHAR
  implicit none

  private
  public :: invLaplace, invHankel, GaussLobatto, TanhSinh, well, formation, solution

  ! Inverse Laplace Transform parameters
  type :: invLaplace
     
     integer :: method = -999
     ! 0 = deHoog, Knight & Stokes
     ! 1 = Stehfest

     ! dehoog-specific parameters
     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha = -999., tol = -999.

     ! 2M+1 Fourier series terms (deHoog)
     ! M terms in Stehfest 
     integer :: M = -999

     ! stehfest-specific salzer summation weights 
     real(EP), allocatable :: V(:)

     ! length of solution vector (2*M+1 or M)
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
     character(81), dimension(9) :: timeDescrip = &
          & ['step on; tpar(1) = on time; tpar(2) not used                                     ',&
          &  'finite pulse; tpar(1:2) = on/off time                                            ',&
          &  'infinitessimal pulse; tpar(1) = pulse location; tpar(2) not used                 ',&
          &  'stairs; tpar(1) = time step (Q increase by integer multiples); tpar(2) = off time',&
          &  'rectified square wave; tpar(1) = 1/2 period of wave; tpar(2) = start time        ',&
          &  'cos(omega*t); tpar(1) = omega; tpar(2) = start time                              ',&
          &  'rectified triangular wave; tpar(1) = 1/4 period of wave; tpar(2) = start time    ',&
          &  'rectified square wave; tpar(1) = 1/2 period of wave; tpar(2) = start time        ',&
          &  'piecewise constant rate (n steps); tpar(1:n)=ti; tpar(n+1)=tfinal; tpar(n+2:)=Q  ']

     ! type of time behavior  (see above)
     integer :: timeType = -999

     ! parameters related to different time behaviors (on, off, etc)
     real(EP), allocatable :: timePar(:)

  end type invLaplace

  ! inverse Hankel transform parameters
  type :: invHankel
     ! zeros of J0 Bessel function
     real(EP), allocatable :: j0z(:) ! locations of zeros of J0 bessel fnc
     integer, allocatable :: sv(:) ! split index  vector

     ! min/max j0 split between infinite/fininte integrals
     integer, dimension(2) :: j0s = [-999, -999]
  end type invHankel

  ! parameters specific to GL quadrature
  type :: GaussLobatto
     integer :: nacc = -999  ! # zeros to accelerate
     integer :: err = -999   ! return value container

     ! abcissa and weights
     real(EP), allocatable :: x(:), w(:)

     ! order of integration
     integer :: ord = -999
  end type GaussLobatto

  type :: vecs
     ! for creating a vector of different-length vectors
     integer, allocatable :: iv(:)
     ! for weights and abcissa at each level
     real(EP), allocatable :: w(:), a(:)
  end type vecs

  ! parameters specific to tanh-sinh quadrature
  type :: TanhSinh
     integer :: k = -999  ! N = 2**k-1 is integration order
     integer :: N = -999
     integer :: R = -999  ! # orders to extrapolate (Rord <= k-2)

     ! these vectors give spacing, order, number and indexing
     ! at each step in the Richardson extrapolation process
     real(EP), allocatable ::  hv(:) ! spacing vector
     integer, allocatable :: kv(:)
     integer, allocatable :: Nv(:)
     type(vecs), allocatable :: Q(:)
  end type TanhSinh

  ! parameters related to well/completion
  type :: well
     real(DP) :: l = -999. ! aquifer top to screen/packer top dist.
     real(DP) :: d = -999. ! aquifer top to screen/packer bottom dist.
     real(DP) :: rw = -999., rc = -999. ! well / casing radii
     real(DP) :: Q = -999. ! volumetric pumping rate

     ! dimensionless parameters
     real(DP) :: rDw = -999. ! dimensionless well radius
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
     ! Mishra & Neuman unsaturated parameters
     ! moisture capacity & hydraulic conductivity sorptive numbers (1/length)
     real(DP) :: ak = -999., ac = -999.
     ! air-entry pressure (<= 0), pressure for saturation (<= 0)
     real(DP) :: psia = 999., psik = 999.
     ! Malama linearization parameter
     real(DP) :: beta = -999.

     real(DP), allocatable :: MoenchAlpha(:), MoenchGamma(:)
     integer :: MoenchAlphaM = -999

     ! computed aquifer parameters
     real(DP) :: sigma = -999.  ! Sy/(Ss*b)
     real(DP) :: alphaD = -999. ! kappa/sigma
     real(DP) :: betaD = -999., usLD = -999.
     real(DP) :: akD = -999., acD = -999., lambdaD = -999.
     real(DP) :: psikD = 999., psiaD = 999., b1 = -999.
     real(DP) :: PsiD = - 999.
  end type formation

  ! parameters related to numerical solution
  type :: solution

     ! characteristic quantities for non-dimensionalizing
     real(DP) :: Lc = -999.  ! length
     real(DP) :: Tc = -999.  ! time
     real(DP) :: Hc = -999.  ! head

     integer :: MNtype = -999  ! switch to pick method of solving Mishra/Neuman2010
     integer :: order = -999 ! order of FD matrix

     ! which unconfined model to use?
     integer :: model = -999
     ! 0 = Theis (confined fully penetrating)
     ! 1 = Hantush (confined partially penetrating)
     ! 2 = Hantush-type with wellbore storage
     ! 3 = Moench 199?
     ! 4 = Malama 2011 fully penetrating
     ! 5 = Malama 2011 partially penetrating
     ! 6 = Mishra/Neuman 2011

     character(15), dimension(0:6) :: modelDescrip = [ &
          & 'Theis          ', 'Hantush        ', 'Hantush w/ stor', &
          & 'Moench         ', 'Malama full pen', 'Malama part pen', &
          & 'Mishra/Neuman  ']

     integer :: quiet = -999  ! verbosity flag (0=quiet, 1=some, 2=lots)
     logical :: dimless = .false.  ! output dimensionless solution?
     logical :: timeseries = .false. ! vector of times, one location?
     logical :: piezometer = .false. ! point observation location?

     ! number of times / locations to be computed or read from file
     integer :: nt = -999, nr = -999, nz = -999

     ! vector of times / locations to compute solution at
     real(DP), allocatable :: t(:), tD(:)
     real(DP), allocatable :: r(:), rD(:)
     real(DP), allocatable :: z(:), zD(:)

     ! which layer (above, in, or below screen) each z point falls in
     integer, allocatable :: zLay(:)

     ! top/bot of monitoring well screen
     real(DP) :: zTop = -999., zBot = -999.
     ! diameter of observation well
     real(DP) :: rwobs = -999., rDwobs = -999.
     ! observation well shape factor
     real(DP) :: sF = -999.
     ! order of quadrature at monitoring well screen
     integer :: zOrd = -999

     character(NUMCHAR) :: outFileName
  end type solution

end module types

