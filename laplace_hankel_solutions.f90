module laplace_hankel_solutions
  implicit none

  private
  public :: lap_hank_soln

contains

  function lap_hank_soln(a,rD,np,nz,w,f,s,lap) result(fp)
    use constants, only : DP, EP
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), intent(in) :: rD
    integer, intent(in) :: np,nz
    type(invLaplace), intent(in) :: lap
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(np,nz) :: fp

    intrinsic :: bessel_j0

    select case(s%model)
    case(0)
       ! Theis solution (no variation in z)
       fp(1:np,1:nz) = spread(2.0/(lap%p(:) + a**2),2,nz)

    case(1)
       ! Hantush solution (confined, partially penetrating)
    case(2)
       ! Boulton model
    case(3)
       ! Neuman 1974 model
    case(4)
       ! Moench 199? model
    case(5)
       ! Mishra/Neuman 2011 model
    case(6)
       ! Malama 2011
    end select

    ! solution always evaluated in test well
    fp(1:np,1:nz) = a*bessel_j0(a*rD)*fp(:,:)*spread(lapTime(lap),2,nz)

  end function lap_hank_soln
    
end module laplace_hankel_solutions




