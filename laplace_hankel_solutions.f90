module laplace_hankel_solutions
  implicit none

  private
  public :: lap_hank_soln

contains

  function lap_hank_soln(a,rD,np,nz,w,f,s,lap) result(fp)
    use constants, only : DP, EP, MAXEXP
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.)  ! outer product operator
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), intent(in) :: rD
    integer, intent(in) :: np,nz
    type(invLaplace), intent(in) :: lap
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(np,nz) :: fp

    complex(EP), allocatable :: eta(:), xi(:), udp(:,:), udf(:,:)

    intrinsic :: bessel_j0

    select case(s%model)
    case(0)
       ! Theis (fully penetrating, confined)
       fp(1:np,1:nz) = theis(a,lap%p,nz)

    case(1)
       ! Hantush (partially penetrating, confined)
       fp(1:np,1:nz) = hantush(a,s%zD,s,lap%p,f,w)

    case(2)
       ! Boulton/Herrera (uconfined, non-physical boundary)
       stop 'ERROR Boulton/Herrera model not implemented yet'
    case(3)
       ! Moench (unconfined w/ wellbore storage)
       stop 'ERROR Moench model not implmented yet'
    case(4)
       ! Malama 2011 fully penetrating model (Neuman 72 when beta=0)
       allocate(eta(np),xi(np),udf(np,nz))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
       udf(1:np,1:nz) = theis(a,lap%p,nz)

       where(spread(abs(eta) < MAXEXP ,2,nz))
          ! naive implementation of formula
          fp(1:np,1:nz) = udf*(1.0_EP - cosh(eta .X. s%zD)/ &
               & spread((1.0_EP + f%beta*eta*xi)*cosh(eta) + xi*sinh(eta),2,nz))
       elsewhere
          ! substitute cosh()&sinh() -> exp() and re-arrange 
          ! only valid when cosh(eta) = 0.5*exp(eta) numerically 
          fp(1:np,1:nz) = udf*(1.0_EP - exp(eta .X. (s%zD - 1.0))/ &
               & spread(1.0_EP + f%betaD*eta*xi + xi,2,nz))
       end where
       deallocate(eta,xi)

    case(5)
       ! Malama 2011 partial penetrating model (Neuman 74 when beta=0)
       allocate(eta(np),xi(np),udp(np,nz+1))

       eta(1:np) = sqrt((lap%p(:) + a**2)/f%kappa)
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
       udp(1:np,1:nz+1) = hantush(a,[s%zD,1.0],s,lap%p,f,w)

       where(spread(abs(eta)<MAXEXP,2,nz))
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & cosh(eta .X. s%zD)/spread((1.0_EP + f%beta*eta*xi)*cosh(eta) + xi*sinh(eta),2,nz)
       elsewhere
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz+1),2,nz)* &
               & exp(eta .X. (s%zD-1.0))/spread(1.0_EP + f%betaD*eta*xi + xi,2,nz)
       end where
       deallocate(eta,xi,udp)

    case(6)
       ! Mishra/Neuman 2011 model
    end select

    fp(1:np,1:nz) = a*bessel_j0(a*rD)*fp(:,:)*spread(lapTime(lap),2,nz)

  end function lap_hank_soln
  
  function theis(a,p,nz) result(udf)
    use constants, only : EP
    real(EP), intent(in) :: a
    complex(EP), dimension(:), intent(in) :: p
    integer, intent(in) :: nz
    complex(EP), dimension(size(p),nz) :: udf

    udf(:,:) = spread(2.0_EP/(p(:) + a**2),dim=2,ncopies=nz)

  end function theis

  function hantush(a,zD,s,p,f,w) result(udp)
    ! implemented in the form given in Malama,Kuhlman & Barrash 2008 
    use constants, only : DP, EP, MAXEXP
    use types, only : well, formation, invLaplace, solution
    use time, only : lapTime
    use utility, only : operator(.X.)  ! outer product operator
    implicit none
    
    real(EP), intent(in) :: a
    real(DP), dimension(:), intent(in) :: zD
    complex(EP), dimension(:), intent(in) :: p
    type(solution), intent(in) :: s
    type(well), intent(in) :: w
    type(formation), intent(in) :: f
    complex(EP), dimension(size(p),size(zd)) :: udp
    complex(EP), dimension(size(p)) :: eta, xi
    complex(EP), dimension(3,size(p),size(zd)) :: g
    integer :: np,nz

    nz = size(zD)
    np = size(p)

    eta(1:np) = sqrt((p(:) + a**2)/f%kappa)
    xi(1:np) = eta(:)*f%alphaD/p(:)

    where(spread(s%zLay(1:nz) == 3,1,np))
       ! above well screen
       g(1,1:np,1:nz) = cosh(eta(1:np) .X. (1.0 - w%dD - zD(1:nz)))
    end where
    
    where(spread(s%zLay(1:nz) == 1,1,np))
       ! below well screen
       g(3,1:np,1:nz) = spread(exp(-eta*(1.0 - w%lD)) - &
            & (sinh(eta*w%dD) + exp(-eta)*sinh(eta*(1.0 - w%lD)))/sinh(eta),2,nz)
    end where

!!$    where(spread(abs(eta(1:np)) < MAXEXP,2,nz))
       g(2,1:np,1:nz) = (spread(sinh(eta*w%dD),2,nz)*cosh(eta .X. zD) + &
                       & spread(sinh(eta*(1.0-w%lD)),2,nz)*cosh(eta .X. (1.0-zD)))/&
                       & spread(sinh(eta),2,nz)
!!$    elsewhere
!!$       g(2,1:np,1:nz) = (exp(eta .X. (w%dD + zD - 1.0)) + &
!!$                       & exp(eta  .X. (w%dD - zd - 1.0)) - &
!!$                       & exp(eta  .X. (zd - w%dD - 1.0)) - &
!!$                       & exp(-eta .X. (w%dD + zd + 1.0)))/2.0_EP + &
!!$                       &(exp(eta  .X. (1.0 - w%lD - zd)) + &
!!$                       & exp(eta  .X. (zd - w%lD - 1.0)) - &
!!$                       & exp(eta  .X. (w%lD - zd - 1.0)) - &
!!$                       & exp(-eta .X. (3.0 - w%lD - zd)))/2.0_EP
!!$    end where

    where(spread(s%zLay(1:nz) == 1,1,np))
       ! below well screen (0 <= zD <= 1-lD)
       udp(1:np,1:nz) =  g(3,:,:)*cosh(eta .X. zD)
    elsewhere
       where(spread(s%zLay(1:nz) == 2,1,np))
          ! next to well screen (1-lD < zD < 1-dD)
          udp(1:np,1:nz) = (1.0_EP - g(2,:,:))
       elsewhere
          ! above well screen (1-lD <= zD <= 1)
          udp(1:np,1:nz) = (g(1,:,:) - g(2,:,:))
       end where
    end where

    udp(1:np,1:nz) = udp(:,:)*theis(a,p,nz)/w%bD
  end function hantush
  
end module laplace_hankel_solutions




