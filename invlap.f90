! this module implements the F.R. de Hoog, J.H. Knight, and A.N. Stokes
! numerical inverse Laplace transform algorithm.
! see "An improved method for numerical inversion of Laplace
!     transforms", SIAM J. Sci. Stat. Comp., 3, 357-366, 1982.

module invlap
  implicit none
  
  private
  public :: deHoog_invlap, deHoog_pvalues
  
  interface deHoog_invlap
     module procedure deHoog_invlap_vect, deHoog_invlap_scalt
  end interface

contains
  
  !! an implementation of the de Hoog et al. method
  !! assumes proper f(p) have been computed for the p
  !! required for the vector of t passed to this function
  !! -- only one log-cycle of time should be passed at once --
  !! (no error checking done in this regard)
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function deHoog_invLap_vect(t,tee,fp,lap) result(ft)
    use constants, only : DP, EP, PI
    use types, only : invLaplace

    real(DP), intent(in) :: tee              ! scaling factor (previously T=2*tmax, but potentially adjustable)
    real(DP), intent(in), dimension(:) :: t   ! vector of times
    type(invLaplace), intent(in) :: lap            ! structure of inputs
    complex(EP), intent(in), dimension(0:2*lap%M) :: fp
    complex(EP), dimension(0:2*lap%M) :: ffpp
    real(EP), dimension(size(t)) :: ft        ! output

    complex(EP), dimension(0:2*lap%M,0:lap%M) :: e
    complex(EP), dimension(0:2*lap%M,1:lap%M) :: q
    complex(EP), dimension(0:2*lap%M) :: d
    complex(EP), dimension(-1:2*lap%M,size(t)) :: A,B
    complex(EP), dimension(size(t)) :: z,brem,rem
    integer ::  r, rq, n, max, nt, M
    real(EP) :: gamma

    M = lap%M
    nt = size(t)

    ! there will be problems if fp(:)==0, or any values are NaN
    if(maxval(abs(fp)) > tiny(1.0_EP)) then

       ffpp = fp
       where(isnan(real(fp)) .or. isnan(aimag(fp)))
          ffpp = (0.0_EP,0.0_EP)
       end where

       ! Re(p) -- this is the de Hoog parameter c
       gamma = lap%alpha - log(lap%tol)/(2.0*tee)

       ! initialize Q-D table 
       e(0:2*M,0) = cmplx(0.0,0.0,EP)
       q(0,1) = ffpp(1)/(ffpp(0)/2.0) ! half first term
       q(1:2*M-1,1) = ffpp(2:2*M)/ffpp(1:2*M-1)

       ! rhombus rule for filling in triangular Q-D table
       do r = 1,M
          ! start with e, column 1, 0:2*M-2
          max = 2*(M-r)
          e(0:max,r) = q(1:max+1,r) - q(0:max,r) + e(1:max+1,r-1)
          if (r /= M) then
             ! start with q, column 2, 0:2*M-3
             rq = r+1
             max = 2*(M-rq)+1
             q(0:max,rq) = q(1:max+1,rq-1) * e(1:max+1,rq-1) / e(0:max,rq-1)
          end if
       end do

       ! build up continued fraction coefficients
       d(0) = ffpp(0)/2.0 ! half first term
       forall(r = 1:M)
          d(2*r-1) = -q(0,r) ! even terms
          d(2*r)   = -e(0,r) ! odd terms
       end forall

       ! seed A and B vectors for recurrence
       A(-1,1:nt) = 0.0
       A(0,1:nt) = d(0)
       B(-1:0,1:nt) = 1.0

       ! base of the power series
       z(1:nt) = exp(cmplx(0.0,1.0,EP)*PI*t(:)/tee)

       ! coefficients of Pade approximation
       ! using recurrence for all but last term
       do n = 1,2*M-1
          A(n,:) = A(n-1,:) + d(n)*A(n-2,:)*z(:)
          B(n,:) = B(n-1,:) + d(n)*B(n-2,:)*z(:)
       end do

       ! "improved remainder" to continued fraction
       brem(1:nt) = (1.0 + (d(2*M-1) - d(2*M))*z(:))/2.0
       rem(1:nt) = -brem*(1.0 - sqrt(1.0 + d(2*M)*z(:)/brem**2))

       ! last term of recurrence using new remainder
       A(2*M,:) = A(2*M-1,:) + rem*A(2*M-2,:)
       B(2*M,:) = B(2*M-1,:) + rem*B(2*M-2,:)

       ! diagonal Pade approximation
       ! F=A/B represents accelerated trapezoid rule
       ft(1:nt) =  exp(gamma*t(:))/tee * real(A(2*M,:)/B(2*M,:))

    else  !! entire f(p) vector is zero
!!$       n = -999
!!$       do r=0,2*M
!!$          if (isnan(abs(ffpp(r)))) then
!!$             n = r
!!$          end if
!!$       end do
!!$       write(*,'(2(A,I0))') 'NaN in f(p), beginning ',n,' out of ',2*M+1
       ft = 0.0
    end if
  end function deHoog_invLap_vect

  function deHoog_invLap_scalt(t,tee,fp,lap) result(ft)
    use constants, only : DP,EP
    use types, only : invLaplace
    real(DP), intent(in) ::  t, tee 
    type(invLaplace), intent(in) :: lap
    complex(EP), intent(in), dimension(0:2*lap%M) :: fp
    real(EP) :: ft ! output
    
    ft = sum(deHoog_invLap_vect([t],tee,fp,lap))
  end function deHoog_invLap_scalt
  
  function deHoog_pvalues(tee,lap) result(p)
    use constants, only : EP, DP, PIEP
    use types, only : invLaplace
    type(invLaplace), intent(in) :: lap
    real(DP), intent(in) :: tee
    complex(EP), dimension(2*lap%M+1) :: p
    real(EP) :: sigma
    integer :: i
    
    ! real portion is constant
    ! TODO: more generally, should the 2.0 in the denominator
    ! TODO: be the constant set in driver.f90?
    sigma = real(lap%alpha,EP) - log(real(lap%tol,EP))/(2.0_EP*tee) 

    forall (i=0:2*lap%M)
       p(i+1) = cmplx(sigma, PIEP*i/tee, EP)
    end forall
    
  end function deHoog_pvalues
end module invlap
