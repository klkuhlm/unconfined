module constants

  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind(p=15,r=300)

  !! extended range internal variables (10 on g95, 10 on gfortran, 16 on ifort)
  integer, parameter :: EP = selected_real_kind(r=3000)

  !! full quad precision (only on gfortran >= 4.6 and ifort)
! integer, parameter :: EP = selected_real_kind(p=33,r=3000)

  !! 3.141592653589793238462643383279503_EP
  real(DP), parameter :: PI =    4.0_DP*atan(1.0_DP) 
  real(EP), parameter :: PIEP =  4.0_EP*atan(1.0_EP) 
  real(EP), parameter :: PIOV2EP = 2.0_EP*atan(1.0_EP)
  real(EP), parameter :: PIOV4EP = atan(1.0_EP)
  complex(EP), parameter :: EYE = (0.0_DP,1.0_EP)
  real(EP), parameter :: E = exp(1.0_EP)
  real(EP), parameter :: SQRT2 = sqrt(2.0_EP)

  !! maximum argument for which sinh(x)-cosh(x) > 0
  ! this is the point where the approximation for sinh()/cosh() -> 0.5*e
  ! for DP  ~18.123, EP ~ 20.72, QP ~ 38.96
  real(EP), parameter :: MAXEXP = -log(epsilon(1.0_EP))/3.0_EP

  !! length of filenames
  integer, parameter :: NUMCHAR = 128

  !! format strings used in output
  character(7) :: RFMT = 'ES14.07'    ! format for general output
  character(9) :: HFMT = 'ES24.15E4'  ! format for results
  character(7) :: SFMT = 'ES09.03'    ! short format for long (+ only) vectors

end module constants
