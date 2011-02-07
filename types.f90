
module types
  use constants, only : DP
  implicit none

  private
  public :: invlt, well, fm, soln

  type :: invlt
     ! Inverse Laplace Transform parameters

     ! abcissa of convergence, LT tolerance
     real(DP) :: alpha = -999., tol = -999.

     ! number of Fourier series terms
     integer :: M = -999
  end type invlt

  type :: well
     ! parameters related to well/completion
     real(DP) :: bD, lD, dD, rDw

  end type well

  type :: fm
     ! parameters related to formation/aquifer
          

     ! dimensionless aquifer parameters
     real(DP) :: alphaD, gamma, kappa, CD
     real(DP), dimension(2) :: beta

  end type fm

  type :: soln
     ! parameters related to numerical solution


  end type soln

end module types

