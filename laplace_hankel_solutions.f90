
module lap_hank_soln
  implicit none

contains

  function laplace_hankel_solution(dum,tD,s,w,f) result(fp)
    use types, only : solution, well, formation

#ifdef INTEL
    use ifport, only : dbesj0 
#endif
    
    implicit none
    
    complex(EP), dimension(s%np) :: soln


    select case()

    ! solution always evaluated in test well
    fp = dum*dbesj0(real(dum*rDw,DP)) * Omega(1:np)

  end function laplace_hankel_solution
  

  function theis() result(fp)

  function unconfined_wellbore_slug(dum,tD,s,w,f) result(fp)
!!$    use shared_data, only :  tsval,bD,dD,lD,alphaD,beta,gamma,rDw,CD,kappa,lap
    use types, only : solution, well, formation
    use constants, only : EP,DP,EONE
    use utilities, only : ccosh, csinh
    use inverse_Laplace_Transform, only : dehoog_pvalues
    use complex_bessel, only : cbesk



    implicit none
    real(EP), intent(in) :: dum  ! scalar integration variable (Hankel parameter)
    real(DP), intent(in) :: tD  ! dimensionless time 
    type(solution) :: s
    type(well) :: w
    type(formation) :: f
    complex(EP), dimension(s%np) :: fp

    complex(EP), dimension(s%np) :: p
    complex(EP), dimension(s%np) :: eta, eps, wD, Omega
    complex(DP), dimension(s%np) :: K1, xiw
    complex(EP), dimension(0:2,s%np) :: delta
    real(EP) :: lDs,dDs
    integer :: i, np,nz,ierr
    logical :: nans

    NaNs = .false.

    np = 2*lap%M+1
    p(1:np) = dehoog_pvalues(2.0*tsval,lap)

    ! pre-compute some intermediate values
    eta(1:np) = sqrt((p(:) + dum**2)/kappa)
    eps(1:np) = p(:)/(eta(:)*alphaD)
    dDs = EONE - real(dD,EP)
    lDs = EONE - real(lD,EP)
    xiw = sqrt(p(:))*rDw
    
    do i=1,np
       call cbesk(z=xiw(i),fnu=1.0D0,kode=1,n=1,cy=K1(i),nz=nz,ierr=ierr)
       if (ierr /= 0 .and. ierr /= 3) then
          write(*,*) 'CBESK error',ierr,' nz=',nz
       end if
    end do

222 continue

    if (maxval(abs(ccosh(eta))) < huge(1.0_EP) .and. (.not. NaNs)) then

       ! when this overflows double precision, switch to approximate form
       delta(0,1:np) = csinh(eta) + eps*ccosh(eta)
       delta(1,1:np) = csinh(eta*dD) + eps*ccosh(eta*dD)
       delta(2,1:np) = csinh(eta*lD) + eps*ccosh(eta*lD)
    
       ! < \hat{ \bar{ w }}_D>
       wD(1:np) = (delta(1,:)*csinh(eta*dDs) + (delta(2,:)-2.0*delta(1,:))*csinh(eta*lds)) / &
            & (bD*eta(:)*delta(0,:))
    else
       ! large argument form (large p or a)
       wD(1:np) = (1.0 - exp(-eta*bD))/(bD*eta)
    end if
    
    ! \hat{ \bar{ \Omega }}
    Omega(1:np) = CD*(EONE - wD(:))/((p(:)+dum**2)*xiw(:)*K1(:))

    if (any(isnan(abs(Omega)))) then
       NaNs = .true.
       goto 222
    end if
    

  end function unconfined_wellbore_slug
end module lap_hank_soln




