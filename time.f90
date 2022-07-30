!
! Copyright (c) 2012-2022 Kristopher L. Kuhlman (klkuhlm at sandia dot gov)
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

module time
  implicit none

  private
  public :: lapTime

contains

  ! ##################################################
  ! return a vector of p for imposing general time behavior
  ! on a function that is a pulse in time
  function lapTime(l) result(mult)
    use constants, only: EP
    use types, only : invLaplace
    implicit none

    type(invLaplace), intent(in) :: l
    complex(EP), dimension(l%np) :: mult

    integer :: n
    real(EP), allocatable :: ti(:), Q(:), W(:), denom(:), y(:)
    real(EP) :: tf

    select case(l%timeType)
    case(1)
       ! step on at time=par1
       mult(1:l%np) = exp(-l%timePar(1)*l%p)/l%p
    case(2)
       ! step on at time=par1, off at time=par2
       mult(1:l%np) = exp(-l%timePar(1)*l%p)/l%p - exp(-l%timePar(2)*l%p)/l%p
    case(3)
       ! instantaneous at t=par1
       mult(1:l%np) = exp(-l%timePar(1)*l%p)
    case(4)
       ! "step test": increasing by integer multiples of Q each
       ! integer multiple of par1 time, off at par2
       mult(1:l%np) = 1.0/(l%p - l%p*exp(-l%timePar(1)*l%p)) * &
                      & (1.0 - exp(-l%timePar(2)*l%p))/l%p
    case(5)
       ! half square wave (only +), period 2*par1
       ! shifted to start at t=par2
       mult(1:l%np) = exp(-l%timePar(2)*l%p)/(l%p + l%p*exp(-l%timePar(1)*l%p))
    case(6)
       ! cosine wave  - cos(at)
       ! shifted to start at t=par2
       mult(1:l%np) = exp(-l%timePar(2)*l%p)*l%p/(l%p**2 + l%timePar(1)**2)
    case(7)
       ! half triangular wave (only +), period 4*par1
       ! shifted to start at t=par2
       mult(1:l%np) = exp(-l%timePar(2)*l%p) / l%p**2 * &
            &(exp(l%timePar(1)*l%p) - exp(l%timePar(1)*l%p))/ &
            &(exp(l%timePar(1)*l%p) + exp(l%timePar(1)*l%p))
    case(8)
       ! full square wave (only +, then -), period 2*par1
       ! shifted to start at t=par2
       mult(1:l%np) = exp(-l%timePar(2)*l%p)* &
            &  (1.0 - exp(-l%timePar(1)*l%p/2.0))/ &
            & ((1.0 + exp(-l%timePar(1)*l%p/2.0))*l%p)
    case(-100:-1)
       !! arbitrary piecewise constant pumping rate with n steps, from ti(1) to tf
       n = -l%timeType
       allocate(ti(n),Q(0:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = l%timePar(1:n)
       tf = l%timePar(n+1)
       Q(0:n) = [0.0_EP, l%timePar(n+2:2*n+1)]

       mult(1:l%np) = (sum(spread(Q(1:n) - Q(0:n-1),2,l%np)*&
            & exp(-spread(ti(1:n),2,l%np)*spread(l%p(1:l%np),1,n)),dim=1) - &
            & sum(Q(1:n) - Q(0:n-1))*exp(-tf*l%p(:)))/l%p(:)

       deallocate(ti,Q)

    case(:-101)
       !! piecewise linear pumping rate with n steps, from ti(1) to tf
       !! no jumps in value (no vertical slopes)
       n = -l%timeType - 100
       allocate(ti(n),W(0:n+1),y(1:n),denom(1:n))

       ! unpack initial times, pumping rates and final time
       ti(1:n) = l%timePar(1:n)
       tf = l%timePar(n+1)
       y(1:n) = l%timePar(n+2:2*n+1)

       ! compute slope between each pair of points
       W(0) = 0.0
       W(n+1) = 0.0
       denom = [ti(2:n),tf] - ti(1:n)
       if (any(abs(denom) < epsilon(1.0))) then
          stop 'no vertical sloped lines in piecewise linear pumping rate'
       end if
       W(1:n) = (y(2:n+1) - y(1:n))/denom ! rise/run

       mult(1:l%np) = (sum(spread(W(1:n) - W(0:n-1),2,l%np)*&
            & exp(spread(-ti,2,l%np)*spread(l%p,1,n)),dim=1) - &
            & sum(W(1:n) - W(0:n-1))*exp(-tf*l%p))/l%p**2

       deallocate(ti,W,y,denom)
    end select

  end function lapTime

end module time
