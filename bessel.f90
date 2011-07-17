SUBROUTINE CJYLV(V,Z,CBJV,CDJV,CBYV,CDYV)
!
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
!
  IMPLICIT DOUBLE PRECISION (A,B,D-H,O-Y)
  IMPLICIT COMPLEX*16 (C,Z)
  integer :: l,k
  DIMENSION CF(12),A(91)
  KM=12
  CALL CJK(KM,A)
  PI=3.141592653589793D0
  DO L=1,0,-1
     V0=V-L
     CWS=SQRT(1.0D0-(Z/V0)*(Z/V0))
     CETA=CWS+LOG(Z/V0/(1.0D0+CWS))
     CT=1.0D0/CWS
     CT2=CT*CT
     DO K=1,KM
        L0=K*(K+1)/2+1
        LF=L0+K
        CF(K)=A(LF)
        DO I=LF-1,L0,-1
           CF(K)=CF(K)*CT2+A(I)
        end DO
        
        CF(K)=CF(K)*CT**K
     end DO
     VR=1.0D0/V0
     CSJ=(1.0D0,0.0D0)
     DO K=1,KM
        CSJ=CSJ+CF(K)*VR**K
     end DO
     
     CBJV=SQRT(CT/(2.0D0*PI*V0))*EXP(V0*CETA)*CSJ
     IF (L.EQ.1) CFJ=CBJV
     CSY=(1.0D0,0.0D0)
     DO K=1,KM
        CSY=CSY+(-1)**K*CF(K)*VR**K
     end DO
        
     CBYV=-SQRT(2.0D0*CT/(PI*V0))*EXP(-V0*CETA)*CSY
     IF (L.EQ.1) CFY=CBYV
  end DO
     
  CDJV=-V/Z*CBJV+CFJ
  CDYV=-V/Z*CBYV+CFY
  RETURN
END SUBROUTINE CJYLV



SUBROUTINE CJK(KM,A)
!
!       ========================================================
!       Purpose: Compute the expansion coefficients for the
!                asymptotic expansion of Bessel functions
!                with large orders
!       Input :  Km   --- Maximum k
!       Output:  A(L) --- Cj(k) where j and k are related to L
!                         by L=j+1+[k*(k+1)]/2; j,k=0,1,...,Km
!       ========================================================
!
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION A(*)
  A(1)=1.0D0
  F0=1.0D0
  G0=1.0D0
  DO K=0,KM-1
     L1=(K+1)*(K+2)/2+1
     L2=(K+1)*(K+2)/2+K+2
     F=(0.5D0*K+0.125D0/(K+1))*F0
     G=-(1.5D0*K+0.625D0/(3.0*(K+1.0D0)))*G0
     A(L1)=F
     A(L2)=G
     F0=F
     G0=G
  end DO
  
  DO K=1,KM-1
     DO J=1,K
        L3=K*(K+1)/2+J+1
        L4=(K+1)*(K+2)/2+J+1
        A(L4)=(J+0.5D0*K+0.125D0/(2.0*J+K+1.0))*A(L3) &
             & -(J+0.5D0*K-1.0+0.625D0/(2.0*J+K+1.0))*A(L3-1)
     end DO
  end DO
  
  RETURN
END SUBROUTINE CJK
