
module laplace_hankel_solution
  implicit none

contains

  function lapl_hank_soln(dum,tD,s,w,f) result(fp)
    use types, only : solution, well, formation
    
    implicit none
    
    complex(EP), dimension(s%np) :: soln


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
    fp = dum*bessel_j0(real(dum*rDw,DP)) * Omega(1:np)

  end function lapl_hank_soln
  

  function theis() result(fp)

  end function theis
  
end module laplace_hankel_solution




