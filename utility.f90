module utility
  implicit none
  private
  public :: logspace, linspace, is_finite, operator(.X.), cosh, sinh

  interface sinh
     module procedure  csinh
  end interface
  
  interface cosh
     module procedure ccosh
  end interface
  
  interface operator(.X.)
     module procedure outerprod_zd, outerprod_dz, outerprod_dd, outerprod_zz
  end interface
  
contains
  
  ! really only needed to make intel compiler work
  elemental function ccosh(z) result(f)
    use constants, only : EP
    complex(EP), intent(in) :: z
    complex(EP) :: f
    real(EP) :: x,y
    intrinsic :: cosh, sinh
    x = real(z)
    y = aimag(z)
    f = cmplx(cosh(x)*cos(y), sinh(x)*sin(y),EP)
  end function ccosh
  
  elemental function csinh(z) result(f)
    use constants, only : EP
    complex(EP), intent(in) :: z
    complex(EP) :: f
    real(EP) :: x,y
    intrinsic :: cosh, sinh
    x = real(z)
    y = aimag(z)
    f = cmplx(sinh(x)*cos(y), cosh(x)*sin(y),EP)
  end function csinh

  function linspace(lo,hi,num) result(v)
    use constants, only : DP
    real(DP), intent(in) :: lo,hi
    integer, intent(in) :: num
    real(DP), dimension(num) :: v
    integer :: i
    real(DP) :: rnum, range, sgn

    rnum = real(num - 1,DP)
    range = abs(hi - lo) 
    sgn = sign(1.0_DP,hi-lo) ! if lo > high, count backwards
    forall (i=0:num-1) 
       v(i+1) = lo + sgn*real(i,DP)*range/rnum
    end forall
  end function linspace

  function logspace(lo,hi,num) result(v)
    use constants, only : DP
    integer, intent(in) :: lo,hi,num
    real(DP), dimension(num) :: v
    v = 10.0_DP**linspace(real(lo,DP),real(hi,DP),num)
  end function logspace

  elemental function is_finite(x) result(pred)
    use constants, only : EP
    complex(EP), intent(in) :: x
    logical :: pred
    pred = .not. (isnan(abs(x)) .or. abs(x) > huge(abs(x)))
  end function is_finite

  pure function outerprod_dd(da,db) result(c)
    use constants, only : DP
    real(DP), intent(in), dimension(:) :: da,db
    real(DP), dimension(size(da),size(db)) :: c
    c = spread(da,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(da))    
  end function outerprod_dd
  
  pure function outerprod_zd(za,db) result(c)
    use constants, only : EP,DP
    complex(EP), intent(in), dimension(:) :: za
    real(DP), intent(in), dimension(:) :: db
    complex(EP), dimension(size(za),size(db)) :: c
    c = spread(za,dim=2,ncopies=size(db))*spread(db,dim=1,ncopies=size(za))
  end function outerprod_zd

  pure function outerprod_dz(da,zb) result(c)
    use constants, only : EP,DP
    real(DP), intent(in), dimension(:) :: da
    complex(EP), intent(in), dimension(:) :: zb
    complex(EP), dimension(size(da),size(zb)) :: c
    c = spread(da,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(da))
  end function outerprod_dz

  pure function outerprod_zz(za,zb) result(c)
    use constants, only : EP
    complex(EP), intent(in), dimension(:) :: za,zb
    complex(EP), dimension(size(za),size(zb)) :: c
    c = spread(za,dim=2,ncopies=size(zb))*spread(zb,dim=1,ncopies=size(za))    
  end function outerprod_zz
  

end module utility


