module laplace_hankel_solutions
  implicit none

  private
  public :: lap_hank_soln

contains

  function lap_hank_soln(a,rD,np,nz,w,f,s,lap) result(fp)
    use constants, only : DP, EP, MAXEXP
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

    complex(EP), allocatable :: eta(:), xi(:)

    intrinsic :: bessel_j0

    select case(s%model)
    case(0)
       ! Theis solution (fully penetrating, confined)
       fp(1:np,1:nz) = theis(a,lap%p,nz)

    case(1)
       ! Hantush solution (confined, partially penetrating)
    case(2)
       ! Boulton model
    case(3)
       ! Neuman 1972 model (fully penetrating, unconfined)
       allocate(eta(np),xi(np))

       eta(1:np) = sqrt(lap%p(:) + a**2)/f%kappa
       xi(1:np) = eta(:)*f%alphaD/lap%p(:)

       ! for large arguments switch to approximation
       where(spread(abs(eta),2,nz) < MAXEXP)
          fp(1:np,1:nz) = theis(a,lap%p,nz)* &
               & (1.0_EP - cosh(spread(eta,2,nz)*spread(s%zD,1,np))/ &
               & spread(cosh(eta(:)) + xi(:)*sinh(eta(:)),2,nz))
       elsewhere
          fp(1:np,1:nz) = theis(a,lap%p,nz)* &
               & (1.0_EP - exp(spread(eta,2,nz)*spread(s%zD - 1.0_EP,1,np))/&
               & spread((1.0_EP + xi(:)),2,nz))
       end where

       deallocate(eta,xi)
    case(4)
!!$       ! Neuman 1974 model (partially penetrating, unconfined)
!!$       allocate(eta(np),xi(np))
!!$
!!$       eta(1:np) = sqrt(lap%p(:) + a**2)/f%kappa
!!$       xi(1:np) = eta(:)*f%alphaD/lap%p(:)
!!$
!!$       ! for large arguments switch to approximation
!!$       where(spread(abs(eta),2,nz) < MAXEXP)
!!$          fp(1:np,1:nz) = spread(2.0_EP/(lap%p(:) + a**2),2,nz)*&
!!$               & (1.0_EP - cosh(spread(eta,2,nz)*spread(s%zD,1,np))/ &
!!$               & spread(cosh(eta(:)) + xi(:)*sinh(eta(:)),2,nz))
!!$       elsewhere
!!$          fp(1:np,1:nz) = spread(2.0_EP/(lap%p(:) + a**2),2,nz)*&
!!$               & (1.0_EP - exp(spread(eta,2,nz)*spread(s%zD - 1.0_EP,1,np))/&
!!$               & spread((1.0_EP + xi(:)),2,nz))
!!$       end where
!!$
!!$       deallocate(eta,xi)
!!$       
!!$       where(spread(s%zLay,1,np) == 1)
!!$          ! below well screen (0 <= zD <= 1-lD)
!!$          udp(1:np,1:nz) = cosh(spread(eta(:),2,nz)*spread(1.0 - w%dD - s%zD,1,np)) - &
!!$               & 
!!$       elsewhere
!!$          where(spread(s%zLay,1,np) == 2)
!!$             ! next to well screen (1-lD < zD < 1-dD)
!!$             udp(1:np,1:nz) = 
!!$          elsewhere
!!$             ! above well screen (1-lD <= zD <= 1)
!!$             udp(1:np,1:nz)
!!$          end where          
!!$       end where
!!$
!!$       udp(1:np,1:nz) = udp*theis(a,p,nz)/w%bD
!!$       
!!$       fp(:,:) = 
!!$
    case(5)
       ! Mishra/Neuman 2011 model
    case(6)
       ! Malama 2011

       

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




