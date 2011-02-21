module laplace_hankel_solutions
  implicit none

  private
  public :: lap_hank_soln

contains

  function lap_hank_soln(a,rD,np,nz,w,frm,s,lap) result(fp)
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
    type(formation), intent(in) :: frm
    complex(EP), dimension(np,nz) :: fp

    integer :: nz1
    real(DP), allocatable :: zd1(:)
    complex(EP), allocatable :: eta(:), xi(:), g(:,:,:), f(:,:), udp(:,:)
    real(EP), allocatable :: abseta(:,:)

    intrinsic :: bessel_j0

    select case(s%model)
    case(0)
       ! Theis (fully penetrating, confined)
       fp(1:np,1:nz) = theis(a,lap%p,nz)
    case(1)
       ! Hantush (partially penetrating, confined)
       stop 'ERROR hantush model not implemented yet'
    case(2)
       ! Boulton/Herrera (uconfined, non-physical boundary)
       stop 'ERROR Boulton/Herrera model not implemented yet'
    case(3)
       ! Moench (unconfined w/ wellbore storage)
       stop 'ERROR Moench model not implmented yet'
    case(4)

       ! Malama 2011 fully penetrating model (Neuman 72 when beta=0)
       allocate(eta(np),xi(np))

       eta(1:np) = sqrt(lap%p(:) + a**2)/frm%kappa
       xi(1:np) = eta(:)*frm%alphaD/lap%p(:)

       where(spread(abs(eta),2,nz) < MAXEXP)
          ! naive implementation of formula
          fp(1:np,1:nz) = theis(a,lap%p,nz)* &
               & (1.0_EP - cosh(eta .X. s%zD)/ &
               & spread(cosh(eta(:)) + xi(:)*sinh(eta(:)),2,nz))
       elsewhere
          ! substitute cosh()&sinh() -> exp() and re-arrange 
          ! only valid when exp(-eta) is numerically = 0.0
          fp(1:np,1:nz) = theis(a,lap%p,nz)* &
               & (1.0_EP - (exp(eta .X. (s%zD - 1.0_DP)) + &
               &           exp(-eta .X. (s%zD + 1.0_DP)))/ &
               & spread(1.0_EP + frm%betaD*eta*xi + xi,2,nz))
       end where
       deallocate(eta,xi)

    case(5)
       ! Malama 2011 partial penetrating model (Neuman 74 when beta=0)
       nz1 = nz+1
       allocate(eta(np),abseta(np,nz1),xi(np),f(3,np),zd1(nz1),udp(np,nz1))
       zd1 = [s%zD,1.0_DP] ! extra zD is for WT boundary (zD=1.0)
       eta(1:np) = sqrt(lap%p(:) + a**2)/frm%kappa
       xi(1:np) = eta(:)*frm%alphaD/lap%p(:)

       f(1,1:np) = sinh(eta(:)*w%dD)
       f(2,1:np) = sinh(eta(:)*(1.0 - w%lD))

       abseta = spread(abs(eta),2,nz1)
       where(abseta(:,1) < MAXEXP)
          ! naive implementation of formula
          f(3,1:np) = exp(-eta(:)*(1.0 - w%lD)) - &
               & (f(1,:) + exp(-eta(:))*f(2,:))/sinh(eta(:))
       elsewhere
          ! substitute cosh()&sinh() -> exp() and re-arrange 
          ! only valid when exp(-eta) is numerically = 0.0
          f(3,1:np) = exp(-eta*w%lD) - exp(eta*(w%dD-1.0_DP)) + &
                    & exp(-eta*(w%dD+1.0_DP)) - exp(-eta*(w%lD+1.0_DP)) + &
                    & exp(-eta*(w%lD+3.0_DP))
       end where
       
       if(any(s%zLay > 1)) then
          ! g not used in layer 1
          allocate(g(2,np,nz1))
          g(1,1:np,1:nz1) = cosh(eta(:) .X. (1.0_DP-w%dD-zD1))

          where(abseta < MAXEXP)
             g(2,1:np,1:nz1) = (spread(f(1,:),2,nz1)*cosh(eta .X. zD1) + &
                              & spread(f(2,:),2,nz1)*cosh(eta .X. (1.0_DP-zD1)))/&
                              & spread(sinh(eta),2,nz1)
          elsewhere
             g(2,1:np,1:nz1) = (exp(eta .X.(w%dD + zD1 - 1.0_DP)) + &
                              & exp(eta .X.(w%dD - zD1 - 1.0_DP)) - &
                              & exp(eta .X.(zD1 - w%dD - 1.0_DP)) - &
                              & exp(-eta.X.(w%dD + zD1 + 1.0_DP)))/2.0_EP + &
                              &(exp(eta .X.(1.0_DP - w%lD - zD1)) + &
                              & exp(eta .X.(zD1 - w%lD - 1.0_DP)) - &
                              & exp(eta .X.(w%lD - zD1 - 1.0_DP)) - &
                              & exp(-eta.X.(3.0_DP - w%lD - zD1)))/2.0_EP
          end where
       end if
       
       where(spread(s%zLay,1,np) == 1)
          ! below well screen (0 <= zD <= 1-lD)
          udp(1:np,1:nz1) =  spread(f(3,:),2,nz1)*cosh(eta .X. zd1) 
       elsewhere
          where(spread(s%zLay,1,np) == 2)
             ! next to well screen (1-lD < zD < 1-dD)
             udp(1:np,1:nz1) = 1.0_EP - g(2,:,:)
          elsewhere
             ! above well screen (1-lD <= zD <= 1)
             udp(1:np,1:nz) = g(1,:,:) - g(2,:,:)
          end where          
       end where
       
       if(any(s%zLay > 1)) then
          deallocate(f,g,zd1)
       else
          deallocate(f,zd1)   
       end if
       
       udp(1:np,1:nz1) = udp(:,:)*theis(a,lap%p,nz1)/w%bD

       where(abseta < MAXEXP)
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz1),2,nz)* &
               & cosh(eta .X. s%zD)/spread(cosh(eta(:)) + xi(:)*sinh(eta(:)),2,nz)
       elsewhere
          fp(1:np,1:nz) = udp(:,1:nz) - spread(udp(:,nz1),2,nz)* &
               & exp(eta .X. (s%zD-1.0))/spread((1.0_EP + xi(:)),2,nz)
       end where

       deallocate(eta,xi,udp)

    case(6)
       ! Mishra/Neuman 2011 model
    end select

    ! solution always evaluated in test well
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

end module laplace_hankel_solutions




