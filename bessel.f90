module complex_bessel

private
public :: cjylv

contains

SUBROUTINE cjylv(v,z,cbjv,cdjv,cbyv,cdyv)

! Code converted using TO_F90 by Alan Miller
! Date: 2011-07-17  Time: 15:26:34

!       ===================================================
!       Purpose: Compute Bessel functions Jv(z) and Yv(z)
!                and their derivatives with a complex
!                argument and a large order
!       Input:   v --- Order of Jv(z) and Yv(z)
!                z --- Complex argument
!       Output:  CBJV --- Jv(z)
!                CDJV --- Jv'(z)
!                CBYV --- Yv(z)
!                CDYV --- Yv'(z)
!       Routine called:
!                CJK to compute the expansion coefficients
!       ===================================================

IMPLICIT DOUBLE PRECISION (a,b,d-h,o-y)
IMPLICIT COMPLEX*16 (c,z)

DOUBLE PRECISION, INTENT(IN)             :: v
COMPLEX, INTENT(IN)                      :: z
COMPLEX, INTENT(OUT)                     :: cbjv
COMPLEX, INTENT(OUT)                     :: cdjv
COMPLEX, INTENT(OUT)                     :: cbyv
COMPLEX, INTENT(OUT)                     :: cdyv
DIMENSION cf(12),a(91)

km=12
CALL cjk(km,a)
pi=3.141592653589793D0
DO  l=1,0,-1
  v0=v-l
  cws=CDSQRT(1.0D0-(z/v0)*(z/v0))
  ceta=cws+CDLOG(z/v0/(1.0D0+cws))
  ct=1.0D0/cws
  ct2=ct*ct
  DO  k=1,km
    l0=k*(k+1)/2+1
    lf=l0+k
    cf(k)=a(lf)
    DO  i=lf-1,l0,-1
      cf(k)=cf(k)*ct2+a(i)
    END DO
    cf(k)=cf(k)*ct**k
  END DO
  vr=1.0D0/v0
  csj=(1.0D0,0.0D0)
  DO  k=1,km
    csj=csj+cf(k)*vr**k
  END DO
  cbjv=CDSQRT(ct/(2.0D0*pi*v0))*CDEXP(v0*ceta)*csj
  IF (l == 1) cfj=cbjv
  csy=(1.0D0,0.0D0)
  DO  k=1,km
    csy=csy+(-1)**k*cf(k)*vr**k
  END DO
  cbyv=-CDSQRT(2.0D0*ct/(pi*v0))*CDEXP(-v0*ceta)*csy
  IF (l == 1) cfy=cbyv
END DO
cdjv=-v/z*cbjv+cfj
cdyv=-v/z*cbyv+cfy
RETURN
END SUBROUTINE cjylv

SUBROUTINE cjk(km,a)

!       ========================================================
!       Purpose: Compute the expansion coefficients for the
!                asymptotic expansion of Bessel functions
!                with large orders
!       Input :  Km   --- Maximum k
!       Output:  A(L) --- Cj(k) where j and k are related to L
!                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
!       ========================================================

IMPLICIT DOUBLE PRECISION (a-h,o-z)
INTEGER, INTENT(IN)                      :: km
DOUBLE PRECISION, INTENT(OUT)            :: a(*)

a(1)=1.0D0
f0=1.0D0
g0=1.0D0
DO  k=0,km-1
  l1=(k+1)*(k+2)/2+1
  l2=(k+1)*(k+2)/2+k+2
  f=(0.5D0*k+0.125D0/(k+1))*f0
  g=-(1.5D0*k+0.625D0/(3.0*(k+1.0D0)))*g0
  a(l1)=f
  a(l2)=g
  f0=f
  g0=g
END DO
DO  k=1,km-1
  DO  j=1,k
    l3=k*(k+1)/2+j+1
    l4=(k+1)*(k+2)/2+j+1
    a(l4)=(j+0.5D0*k+0.125D0/(2.0*j+k+1.0))*a(l3)  &
        -(j+0.5D0*k-1.0+0.625D0/(2.0*j+k+1.0))*a(l3-1)
  END DO
END DO
RETURN
END SUBROUTINE cjk


end module complex_bessel
