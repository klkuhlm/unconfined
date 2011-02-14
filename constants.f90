module constants

  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind (p=15,r=300)

  !! extended range internal variables (10 on g95, 10 on gfortran, 16 on ifort)
  integer, parameter :: EP = selected_real_kind(r=3000)

  !! full quad precision (only on gfortran >= 4.6 and ifort)
  integer, parameter :: QP = selected_real_kind(p=33,r=3000)

  !! 3.141592653589793238462643383279503_EP
  real(DP), parameter :: PI =    4.0_DP*atan(1.0_DP) 
  real(EP), parameter :: PIEP =  4.0_EP*atan(1.0_EP) 
  real(DP), parameter :: PIOV2EP = 2.0_EP*atan(1.0_EP)

  !! length of filenames
  integer, parameter :: NUMCHAR = 128

end module constants
