module utility
  implicit none
  private
  public :: logspace, linspace, is_finite

contains
  function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    integer :: i
    real(DP) :: rnum, range, sgn

    rnum = real(num - 1,DP)
    range = abs(hi - lo) 
    sgn = sign(1.0_DP,hi-lo) ! if lo > high, count backwards
    
    !$OMP PARALLEL WORKSHARE
    forall (i=0:num-1) 
       v(i+1) = lo + sgn*real(i,DP)*range/rnum
    end forall
    !$OMP END PARALLEL WORKSHARE

  end function linspace

  function logspace(lo,hi,num) result(v)
    use constants, only : DP
    integer, intent(in) :: lo,hi,num
    real(DP), dimension(num) :: v
    v = 10.0_DP**linspace(real(lo,DP),real(hi,DP),num)
  end function logspace

  elemental function is_finite(x) result(pred)
    use constants, only : EP
    complex(EP), intent(in) :: x
    logical :: pred
    pred = .not. (isnan(abs(x)) .or. abs(x) > huge(abs(x)))
  end function is_finite

end module utility


