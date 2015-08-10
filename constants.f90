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

module constants

  public

  ! real with range 300 orders of mag, 15 sig figs (8 on gfortran, g95 & ifort)
  integer, parameter :: DP = selected_real_kind(p=15,r=300)

  ! changing the definition of EP ("extended precision") to one of the
  ! following 3 options will allow the code to run faster or be more
  ! precise.

  ! 1) Setting EP=8 or DP (double precision) is the fastest.
  integer, parameter :: EP = DP

  ! 2) Setting EP=10 has a bit more precision but quad-precision range
  ! (implemented in hardward on intel/amd chipsets with 80-bit
  ! registers).

  !! extended range internal variables (10 on g95, 10 on gfortran, 16 on ifort)
  !! integer, parameter :: EP = selected_real_kind(r=3000)
  
  ! 3) Setting EP=16 is full quad-precision range and precision, but
  ! is very slow (implemented in software)

  !! full quad precision (only on gfortran >= 4.6 and ifort)
  integer, parameter :: QP = selected_real_kind(p=33,r=3000)

  !! 3.141592653589793238462643383279503_EP
  real(DP), parameter :: PI =    4.0_DP*atan(1.0_DP)
  real(DP), parameter :: LN2 = log(2.0_DP)

  real(EP), parameter :: PIEP =  4.0_EP*atan(1.0_EP)
  real(EP), parameter :: INVPIEP = 0.25_EP/atan(1.0_EP)
  real(EP), parameter :: TWOPIEP = 8.0_EP*atan(1.0_EP)
  real(EP), parameter :: PIOV2EP = 2.0_EP*atan(1.0_EP)
  real(EP), parameter :: PIOV4EP = atan(1.0_EP)
  complex(EP), parameter :: EYE = (0.0_DP,1.0_EP)
  real(EP), parameter :: E = exp(1.0_EP)
  real(EP), parameter :: SQRT2 = sqrt(2.0_EP)

  !! maximum argument for which sinh(x)-cosh(x) > 0
  ! this is the point where the approximation for sinh(x) or cosh(x) -> 0.5*exp(x)
  ! for DP ~ 18.123, EP ~ 20.72, QP ~ 38.96
  real(EP), parameter :: MAXEXP = -log(epsilon(1.0_EP))/3.0_EP

  !! length of filenames
  integer, parameter :: NUMCHAR = 128

  !! format strings used in output
  character(9) :: RFMT = 'ES14.07E2'    ! format for general output
  character(9) :: HFMT = 'ES24.15E4'    ! format for results
  character(9) :: SFMT = 'ES09.03E2'    ! short format for long (+ only) vectors

  ! NB: changing these formats will mess with column alignment in
  ! output (but this is just an astetic)

end module constants
