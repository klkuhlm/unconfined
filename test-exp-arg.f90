program testexp
  implicit none

  ! real with range 300 orders of mag, 15 sig figs (8 on both g95 & ifort)
  integer, parameter :: DP = selected_real_kind(p=15,r=300)

  !! extended range internal variables (10 on g95, 10 on gfortran, 16 on ifort)
  integer, parameter :: EP = selected_real_kind(r=3000)

  !! full quad precision (only on gfortran >= 4.6 and ifort)
  integer, parameter :: QP = selected_real_kind(p=33,r=3000)

  real(DP), parameter :: LOG2DP = log10(2.0_DP)
  real(EP), parameter :: LOG2EP = log10(2.0_EP)
  real(QP), parameter :: LOG2QP = log10(2.0_QP)

  ! complex should behave like reals of same kind
  real(DP) :: x,y!,z
  real(EP) :: i,j!,k
  real(QP) :: a,b!,c

  integer :: ii
  character(38) :: fmt

  ! look at the numerical behavior of cosh/sinh/exp for large
  ! arguments.  Need to figure out the optimum/appropirate
  ! changeovers for the naive and large-argument forms.

  ! Looking for two points.  First:
  ! the point where cosh(x)-sinh(x) = 0.0 numerically.
  ! Secondly, where exp(-x) = 0.0 numerically.  These are NOT the
  ! same, as I initially assumed they were.

  ! the first case occurs for a rather small x, due to the 
  ! subtraction of like quantities and severe cancellation.
  ! the second case is much larger.

  ! when x is beyond the size of the second case, the large-argument
  ! form is the most accurate form available.  Because
  ! cosh(x) = 0.5*(exp(x) + exp(-x)) and
  ! sinh(x) = 0.5*(exp(x) - exp(-x)), in this region it is
  ! true that cosh(x) = sinh(x) = 0.5*exp(x).  This should be handled
  ! properly by the library implementation of cosh and sinh,
  ! but can often lead to significant simplification of the
  ! laplace-hankel space argument.

  ! in the intermediate region (smaller than the second case, but
  ! larger than the first case), re-arrangement leads to the
  ! obvious tautology that cosh(x)-sinh(x)=exp(-x), but numerically
  ! the cancellation makes this zero for much smaller argument.
  ! This is only important for cases where there are subtraction of
  ! sinh and cosh of identical or similar argument.

  ! The smaller number, reflected in the underflow of the exponential
  ! for very large argument.  Determine the cutoff where this should
  ! happen for a given real/complex datatype (single-,double-,
  ! quad-precision).

  fmt = '(I3,1X,ES22.15,2(1X,I3),2(1X,ES22.15))   '
  x = 1.0; i = 1.0; a = 1.0
  print *, '-log(epsilon(1.0_DP))',-log(epsilon(x)) ! ~36.04
  print *, 'digits(1.0_DP)       ',digits(x) ! 2^53
  print *, 'precision(1.0_DP)    ',precision(x),(digits(x)-1)*log10(2.0_DP) ! 10^15
  print *, '{min,max}exponent(1.0_DP)',minexponent(x),maxexponent(x)
  print *, 'epsilon(1.0_DP)',epsilon(x),2.0_DP**(1-digits(x)),2.0_DP**(-digits(x))
  print *, '1-eps',1.0_DP-2.0_DP**(1-digits(x)),1.0_DP-2.0_DP**(-digits(x))
  print *, '10**(1-precision(x))',10.0_DP**(1-(digits(x)-1)*log10(2.0_DP)) ! _not_ same as epsilon
  print *, 'log10(2.0)',log10(2.0_DP)
  print *, 'huge(x)',huge(x),1.0_EP-2.0_DP**(-digits(x))
  do ii = 0,precision(x)+2
     y = 1.0_DP - (1.0_DP + 10.0_DP**(-ii))
     write(*,fmt) ii,y,-exponent(y),-int(exponent(y)*LOG2EP),fraction(1.0_DP),fraction((1.0_DP + 10.0_DP**(-ii)))
  end do
  
  fmt = '(I3,1X,ES25.18,2(1X,I3),2(1X,ES25.18))               '
  print *, ' '
  print *, '-log(epsilon(1.0_EP))',-log(epsilon(i)) ! ~43.67
  print *, 'digits(1.0_EP)       ',digits(i) ! 2^64
  print *, 'precision(1.0_EP)    ',precision(i),(digits(i)-1)*log10(2.0_EP) ! 10^18
  print *, '{min,max}exponent(1.0_EP)',minexponent(i),maxexponent(i)
  print *, 'epsilon(1.0_EP)',epsilon(i),2.0_EP**(1-digits(i)),2.0_EP**(-digits(i))
  print *, '1-eps',1.0_EP-2.0_EP**(1-digits(i)),1.0_EP-2.0_EP**(-digits(i))
  print *, '10**(1-precision(x))',10.0_EP**(1-precision(i))
  print *, 'log10(2.0)',log10(2.0_EP)
  print *, 'huge(x)',huge(i),1.0_EP-2.0_EP**(-digits(i))
  do ii=0,precision(i)+2
     j = 1.0_EP - (1.0_EP + 10.0_EP**(-ii))
     write(*,fmt) ii,j,-exponent(j),-int(exponent(j)*LOG2EP),fraction(1.0_EP),fraction((1.0_EP + 10.0_EP**(-ii)))
  end do

  fmt = '(I3,1X,ES40.33,2(1X,I3),2(1X,ES40.33))               '
  print *, ' '
  print *, '-log(epsilon(1.0_QP))',-log(epsilon(a)) ! ~77.63
  print *, 'digits(1.0_QP)       ',digits(a) ! 2^113
  print *, 'precision(1.0_QP)    ',precision(a),(digits(a)-1)*log10(2.0_QP) ! 10^33
  print *, '{min,max}exponent(1.0_QP)',minexponent(a),maxexponent(a)
  print *, 'epsilon(1.0_QP)',epsilon(a),2.0_QP**(1-digits(a)),2.0_QP**(-digits(a))
  print *, '1-eps',1.0_QP-2.0_QP**(1-digits(a)),1.0_QP-2.0_QP**(-digits(a))
  print *, '10**(1-precision(x))',10.0_QP**(1-precision(a))
  print *, 'log10(2.0)',log10(2.0_QP)
  print *, 'huge(x)',huge(a),1.0_QP-2.0_QP**(-digits(a))
  do ii=0,precision(a)+2
     b = 1.0_QP - (1.0_QP + 10.0_QP**(-ii))
     write(*,fmt) ii,b,-exponent(b),-int(exponent(b)*LOG2QP),fraction(1.0_QP),fraction((1.0_QP + 10.0_QP**(-ii)))
  end do

end program testexp
