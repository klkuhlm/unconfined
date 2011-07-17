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
    use constants, only : EP,PIEP
    implicit none

    REAL(EP), INTENT(IN)  :: v
    COMPLEX(EP), INTENT(IN)   :: z
    COMPLEX(EP), INTENT(OUT)  :: cbjv
    COMPLEX(EP), INTENT(OUT)  :: cdjv
    COMPLEX(EP), INTENT(OUT)  :: cbyv
    COMPLEX(EP), INTENT(OUT)  :: cdyv

    real(EP), dimension(91) :: a
    complex(EP), dimension(12) :: cf

    integer :: i,l,k,km,l0,lf
    complex(EP) :: cws,ceta,ct,ct2,csj,csy,cfj,cfy
    real(EP) :: v0,vr

    km=12
    CALL cjk(km,a)
    DO l=1,0,-1
       v0 = v-l
       cws = SQRT(1.0_EP-(z/v0)*(z/v0))
       ceta = cws+LOG(z/v0/(1.0_EP+cws))
       ct = 1.0_EP/cws
       ct2 = ct*ct
       DO k=1,km
          l0 = k*(k+1)/2+1
          lf = l0+k
          cf(k) = a(lf)
          DO i=lf-1,l0,-1
             cf(k) = cf(k)*ct2+a(i)
          END DO
          cf(k) = cf(k)*ct**k
       END DO
       vr = 1.0_EP/v0
       csj = (1.0_EP,0.0_EP)
       DO k=1,km
          csj = csj+cf(k)*vr**k
       END DO
       cbjv = SQRT(ct/(2.0_EP*PIEP*v0))*EXP(v0*ceta)*csj
       IF (l == 1) then
          cfj = cbjv
       end IF
       csy = (1.0_EP,0.0_EP)
       DO k=1,km
          csy = csy+(-1)**k*cf(k)*vr**k
       END DO
       cbyv = -SQRT(2.0_EP*ct/(PIEP*v0))*EXP(-v0*ceta)*csy
       IF (l == 1) then
          cfy=cbyv
       end IF
    END DO
    cdjv = -v/z*cbjv+cfj
    cdyv = -v/z*cbyv+cfy
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
    use constants, only : EP
    implicit none
    INTEGER, INTENT(IN) :: km
    real(EP), INTENT(OUT) :: a(*)

    integer :: k,j,l1,l2,l3,l4
    real(EP) :: f0,g0,f,g

    a(1) = 1.0_EP
    f0 = 1.0_EP
    g0 = 1.0_EP
    DO k=0,km-1
       l1 = (k+1)*(k+2)/2+1
       l2 = (k+1)*(k+2)/2+k+2
       f = (0.5_EP*k+0.125_EP/(k+1))*f0
       g = -(1.5_EP*k+0.625_EP/(3.0_EP*(k+1.0_EP)))*g0
       a(l1) = f
       a(l2) = g
       f0 = f
       g0 = g
    END DO
    DO k=1,km-1
       DO j=1,k
          l3 = k*(k+1)/2+j+1
          l4 = (k+1)*(k+2)/2+j+1
          a(l4) = (j+0.5_EP*k+0.125_EP/(2.0_EP*j+k+1.0_EP))*a(l3)  &
               -(j+0.5_EP*k-1.0_EP+0.625_EP/(2.0_EP*j+k+1.0_EP))*a(l3-1)
       END DO
    END DO
  END SUBROUTINE cjk
end module complex_bessel
