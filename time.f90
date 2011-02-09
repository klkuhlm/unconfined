module time_mod
  implicit none

  private 
  public :: time

contains

  ! ##################################################
  ! general time behavior function, for all elements
   pure function time(p,t) result(mult)
    use constants, only: EP
    use types, only : solution
    implicit none

    type(solution), intent(in) :: s
    complex(EP), dimension(:), intent(in) :: p
    complex(EP), dimension(size(p,1)) :: mult

    integer :: n
    real(EP), allocatable :: ti(:), Q(:)
    real(EP) :: tf
    
    select case(s%timeType)
    case(1)
       ! step on at time=par1
       mult(1:s%np) = exp(-s%timePar(1)*p)/p
    case(2)
       ! step on at time=par1, off at time=par2
       mult(1:s%np) = exp(-s%timePar(1)*p)/p - exp(-s%timePar(2)*p)/p
    case(3)
       ! instantaneous at t=par1
       mult(1:s%np) = exp(-s%timePar(1)*p)
    case(4)  
       ! "step test": increasing by integer multiples of Q each 
       ! integer multiple of par1 time, off at par2
       mult(1:s%np) = 1.0/(p - p*exp(-s%timePar(1)*p)) * &
                      & (1.0 - exp(-s%timePar(2)*p))/p
    case(5)
       ! half square wave (only +), period 2*par1
       ! shifted to start at t=par2
       mult(1:s%np) = exp(-s%timePar(2)*p)/(p + p*exp(-s%timePar(1)*p))
    case(6)
       ! cosine wave  - cos(at)
       ! shifted to start at t=par2
       mult(1:s%np) = exp(-s%timePar(2)*p)*p/(p**2 + s%timePar(1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       ! shifted to start at t=par2
       mult(1:s%np) = exp(-s%timePar(2)*p) / p**2 * &
            &(exp(s%timePar(1)*p) - exp(s%timePar(1)*p))/ &
            &(exp(s%timePar(1)*p) + exp(s%timePar(1)*p))
    case(8)
       ! full square wave (only +, then -), period 2*par1
       ! shifted to start at t=par2
       mult(1:s%np) = exp(-s%timePar(2)*p)* &
            &  (1.0 - exp(-s%timePar(1)*p/2.0))/ &
            & ((1.0 + exp(-s%timePar(1)*p/2.0))*p)
    case(:-1)
       !! arbitrary piecewise constant pumping rate with n steps, from ti(1) to tf
       n = -s%timeType
       allocate(ti(n),Q(0:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = s%timePar(1:n)
       tf = s%timePar(n+1)
       Q(0:n) = [0.0, s%timePar(n+2:2*n+1)]
       
       mult(1:s%np) = (sum(spread(Q(1:n) - Q(0:n-1),2,s%np)*&
            & exp(-spread(ti(1:n),2,s%np)*spread(p(1:s%np),1,n)),dim=1) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*p(:)))/p(:)

       deallocate(ti,Q)
    end select
    deallocate(par)

  end function Time

end module time_mod
