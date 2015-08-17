program drive_arb

  ! compute bessel function using Arb

  use, intrinsic :: iso_c_binding, only : C_FLOAT128_COMPLEX, C_FLOAT128
  
  implicit none
  interface
     function arb_J(nu,z) bind(c,name="arb_J") result(J)
       use, intrinsic :: iso_c_binding, only : C_FLOAT128_COMPLEX, C_FLOAT128
       complex(C_FLOAT128_COMPLEX), intent(in), value :: z
       real(C_FLOAT128), intent(in), value :: nu
       complex(C_FLOAT128_COMPLEX) :: J
     end function arb_J
  end interface

  interface
     function arb_Y(nu,z) bind(c,name="arb_Y") result(Y)
       use, intrinsic :: iso_c_binding, only : C_FLOAT128_COMPLEX, C_FLOAT128
       complex(C_FLOAT128_COMPLEX), intent(in), value :: z
       real(C_FLOAT128), intent(in), value :: nu
       complex(C_FLOAT128_COMPLEX) :: Y
     end function arb_Y
  end interface
  
  complex(C_FLOAT128_COMPLEX) :: z
  real(C_FLOAT128_COMPLEX), parameter :: PIOV4 = atan(1.0_C_FLOAT128_COMPLEX)
  real(C_FLOAT128) :: nu
  integer :: i

  do i=1,8
     z = cmplx(i*PIOV4, 2*i*PIOV4, kind=C_FLOAT128_COMPLEX)
     nu = (i+2)*PIOV4
     print *, i,z,nu,arb_J(nu,z),arb_Y(nu,z)
  end do

end program drive_arb
