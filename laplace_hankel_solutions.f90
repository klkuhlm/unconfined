
module laplace_hankel_solution
  implicit none

contains

  function lapl_hank_soln(a,tD,rD,np,nz,w,f,lap) result(fp)
    use constants, only : DP, EP
    use types, only : well, formation, invLaplace
    use time_mod, only : lapTime
    
    implicit none
    
    real(DP), intent(in) :: a
    real(DP), intent(in) :: tD,rD
    integer, intent(in) :: np,nz
    type(invLaplace), intent(in) :: lap
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(np,nz) :: fp

    select case(s%model)
    case(0)
       ! Theis solution
       
       

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
    fp = a*bessel_j0(a*rD)*fp*spread(lapTime(lap),2,nz)

  end function lapl_hank_soln
    
end module laplace_hankel_solution




