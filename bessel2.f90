module complex_bessel2

private
public :: cjyva

contains
SUBROUTINE cjyva(v,z,vm,cbj,cdj,cby,cdy)

  !       ===========================================================
  !       Purpose: Compute Bessel functions Jv(z), Yv(z) and their
  !                derivatives for a complex argument
  !       Input :  z --- Complex argument
  !                v --- Order of Jv(z) and Yv(z)
  !                      ( v = n+v0, n = 0,1,2,..., 0 ≤ v0 < 1 )
  !       Output:  CBJ(n) --- Jn+v0(z)
  !                CDJ(n) --- Jn+v0'(z)
  !                CBY(n) --- Yn+v0(z)
  !                CDY(n) --- Yn+v0'(z)
  !                VM --- Highest order computed
  !       Routines called:
  !            (1) GAMMA2 for computing the gamma function
  !            (2) MSTA1 and MSTA2 for computing the starting
  !                point for backward recurrence
  !       ===========================================================

  use constants, only : EP, PI

  IMPLICIT real(EP) (a,b,g,o-y)
  IMPLICIT complex(EP) (c,z)

  real(EP), INTENT(IN) :: v
  complex(EP), INTENT(IN) :: z
  real(EP), INTENT(OUT) :: vm
  COMPLEX(EP), INTENT(OUT) :: cbj(0:int(v)+1)
  COMPLEX(EP), INTENT(OUT) :: cdj(0:int(v)+1)
  COMPLEX(EP), INTENT(OUT) :: cby(0:int(v)+1)
  COMPLEX(EP), INTENT(OUT) :: cdy(0:int(v)+1)

  intrinsic :: gamma

  rp2=.63661977236758_EP
  ci=(0.0_EP,1.0_EP)
  a0=ABS(z)
  z1=z
  z2=z**2
  n=INT(v)
  v0=v-n
  pv0=PI*v0
  pv1=PI*(1.0_EP+v0)
  ! Treat z=0 as a special case
  IF (a0 < 1.0D-100) THEN
     cbj(0:n)=(0.0_EP,0.0_EP)
     cdj(0:n)=(0.0_EP,0.0_EP)
     cby(0:n)=-(1.0D-300,0.0_EP)
     cdy(0:n)=(1.0D-300,0.0_EP)
     IF (v0 == 0.0) THEN
        cbj(0)=(1.0_EP,0.0_EP)
        cdj(1)=(0.5_EP,0.0_EP)
     ELSE
        cdj(0)=(1.0D-300,0.0_EP)
     END IF
     vm=v
     RETURN
  END IF
  
  lb0=0
  IF (real(z) < 0.0) z1=-z
  ! Calculate J_{v_0} and J_{v_0+1} using eqn (5.5.1) for |z|<=12
  IF (a0 <= 12.0) THEN
     DO  l=0,1
        vl=v0+l
        cjvl=(1.0_EP,0.0_EP)
        cr=(1.0_EP,0.0_EP)
        loopa: DO  k=1,40
           cr=-0.25_EP*cr*z2/(k*(k+vl))
           cjvl=cjvl+cr
           IF (ABS(cr) < ABS(cjvl)*1.0D-15) then 
              EXIT loopa
           end IF
        END DO loopa
        vg=1.0_EP+vl
        ga = gamma(vg)
        ca=(0.5_EP*z1)**vl/ga
        IF (l == 0) cjv0=cjvl*ca
        IF (l == 1) cjv1=cjvl*ca
     END DO
  ELSE
     ! otherwise, use eqns (5.2.5) and (5.2.6)
     k0=11
     IF (a0 >= 35.0) k0=10
     IF (a0 >= 50.0) k0=8
     DO  j=0,1
        vv=4.0_EP*(j+v0)*(j+v0)
        cpz=(1.0_EP,0.0_EP)
        crp=(1.0_EP,0.0_EP)
        DO  k=1,k0
           crp=-0.0078125_EP*crp*(vv-(4.0*k-3.0)**2.0)*(vv-  &
                (4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
           cpz=cpz+crp
        END DO
        cqz=(1.0_EP,0.0_EP)
        crq=(1.0_EP,0.0_EP)
        DO  k=1,k0
           crq=-0.0078125_EP*crq*(vv-(4.0*k-1.0)**2.0)*(vv-  &
                (4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
           cqz=cqz+crq
        END DO
        cqz=0.125_EP*(vv-1.0)*cqz/z1
        zk=z1-(0.5_EP*(j+v0)+0.25_EP)*pi
        ca0=SQRT(rp2/z1)
        cck=COS(zk)
        csk=SIN(zk)
        ! Note that Y_{v_0} and Y_{v_0+1} are also calculated here for |z|<=12
        IF (j == 0) THEN
           cjv0=ca0*(cpz*cck-cqz*csk)
           cyv0=ca0*(cpz*csk+cqz*cck)
        ELSE IF (j == 1) THEN
           cjv1=ca0*(cpz*cck-cqz*csk)
           cyv1=ca0*(cpz*csk+cqz*cck)
        END IF
     END DO
  END IF
  IF (a0 <= 12.0) THEN
     IF (v0 /= 0.0) THEN
        ! Calculate J_{-v_0} and J_{-(v_0+1)} using eqn (5.5.1) for |z|<=12
        DO  l=0,1
           vl=v0+l
           cjvl=(1.0_EP,0.0_EP)
           cr=(1.0_EP,0.0_EP)
           loopb: DO  k=1,40
              cr=-0.25_EP*cr*z2/(k*(k-vl))
              cjvl=cjvl+cr
              IF (ABS(cr) < ABS(cjvl)*1.0D-15) then
                 EXIT loopb
              end IF
           END DO loopb
           vg=1.0_EP-vl
           gb = gamma(vg)
           cb=(2.0_EP/z1)**vl/gb
           IF (l == 0) cju0=cjvl*cb
           IF (l == 1) cju1=cjvl*cb
        END DO
        ! Calculate Y_{v_0} and Y_{v0+1} using eqn (5.1.4)
        cyv0=(cjv0*COS(pv0)-cju0)/SIN(pv0)
        cyv1=(cjv1*COS(pv1)-cju1)/SIN(pv1)
     ELSE
        cec=LOG(z1/2.0_EP)+.5772156649015329_EP
        ! Calculate Y_0 and Y_1 using eqns (5.2.3) and (5.2.4)
        cs0=(0.0_EP,0.0_EP)
        w0=0.0_EP
        cr0=(1.0_EP,0.0_EP)
        DO  k=1,30
           w0=w0+1.0_EP/k
           cr0=-0.25_EP*cr0/(k*k)*z2
           cs0=cs0+cr0*w0
        END DO
        cyv0=rp2*(cec*cjv0-cs0)
        cs1=(1.0_EP,0.0_EP)
        w1=0.0_EP
        cr1=(1.0_EP,0.0_EP)
        DO  k=1,30
           w1=w1+1.0_EP/k
           cr1=-0.25_EP*cr1/(k*(k+1))*z2
           cs1=cs1+cr1*(2.0_EP*w1+1.0_EP/(k+1.0_EP))
        END DO
        cyv1=rp2*(cec*cjv1-1.0_EP/z1-0.25_EP*z1*cs1)
     END IF
  END IF
  IF (DBLE(z) < 0.0_EP) THEN
     ! Apply eqns (5.4.1) and (5.4.2)
     cfac0=EXP(pv0*ci)
     cfac1=EXP(pv1*ci)
     IF (AIMAG(z) < 0.0_EP) THEN
        cyv0=cfac0*cyv0-2.0_EP*ci*COS(pv0)*cjv0
        cyv1=cfac1*cyv1-2.0_EP*ci*COS(pv1)*cjv1
        cjv0=cjv0/cfac0
        cjv1=cjv1/cfac1
     ELSE IF (AIMAG(z) > 0.0_EP) THEN
        cyv0=cyv0/cfac0+2.0_EP*ci*COS(pv0)*cjv0
        cyv1=cyv1/cfac1+2.0_EP*ci*COS(pv1)*cjv1
        cjv0=cfac0*cjv0
        cjv1=cfac1*cjv1
     END IF
  END IF
  cbj(0)=cjv0
  cbj(1)=cjv1
  IF (n >= 2.AND.n <= INT(0.25*a0)) THEN
     ! Calculate J_v using forward recurrence for N<|z|/4
     cf0=cjv0
     cf1=cjv1
     DO  k=2,n
        cf=2.0_EP*(k+v0-1.0_EP)/z*cf1-cf0
        cbj(k)=cf
        cf0=cf1
        cf1=cf
     END DO
  ELSE IF (n >= 2) THEN
     ! otherwise use backwards recurrence
     m=msta1(a0,2000) ! increased from 200 -> 1000
     IF (m < n) THEN
        n=m
     ELSE
        if (EP == 10) then
           m=msta2(a0,n,17) ! increased from 15 -> 17 (kind = 10)
        else
           m=msta2(a0,n,30) ! increased from 15 -> 30 (kind = 16)
        end if
     END IF
     cf2=(0.0_EP,0.0_EP)
     cf1=(1.0D-100,0.0_EP)
     DO  k=m,0,-1
        cf=2.0_EP*(v0+k+1.0_EP)/z*cf1-cf2
        IF (k <= n) cbj(k)=cf
        cf2=cf1
        cf1=cf
     END DO
     IF (ABS(cjv0) > ABS(cjv1)) cs=cjv0/cf
     IF (ABS(cjv0) <= ABS(cjv1)) cs=cjv1/cf2
     DO  k=0,n
        cbj(k)=cs*cbj(k)
     END DO
  END IF
  cdj(0)=v0/z*cbj(0)-cbj(1)
  DO  k=1,n
     ! calculate the derivative of J_v
     cdj(k)=-(k+v0)/z*cbj(k)+cbj(k-1)
  END DO
  cby(0)=cyv0
  cby(1)=cyv1
  ! NOTE: the following part calculates Y_v using mixed recurrence.  It
  ! can be replaced by a simpler algorithm based on eqns (5.5.11) and (5.5.12)
  ya0=ABS(cyv0)
  lb=0
  ! determine the turning point B approximately using forward recurrence
  cg0=cyv0
  cg1=cyv1
  DO  k=2,n
     cyk=2.0_EP*(v0+k-1.0_EP)/z*cg1-cg0
     IF (ABS(cyk) > 1.0D+290) CYCLE
     yak=ABS(cyk)
     ya1=ABS(cg0)
     IF (yak < ya0.AND.yak < ya1) lb=k
     cby(k)=cyk
     cg0=cg1
     cg1=cyk
  END DO
  IF (lb <= 4.OR.AIMAG(z) == 0.0_EP) GO TO 125
95 IF (lb == lb0) GO TO 125
  ch2=(1.0_EP,0.0_EP)
  ch1=(0.0_EP,0.0_EP)
  lb0=lb
  DO  k=lb,1,-1
     ch0=2.0_EP*(k+v0)/z*ch1-ch2
     ch2=ch1
     ch1=ch0
  END DO
  ! Calculate P_{12} and P_{22}
  cp12=ch0
  cp22=ch2
  ch2=(0.0_EP,0.0_EP)
  ch1=(1.0_EP,0.0_EP)
  DO  k=lb,1,-1
     ch0=2.0_EP*(k+v0)/z*ch1-ch2
     ch2=ch1
     ch1=ch0
  END DO
  ! Calculate P_{11} and P_{21}
  cp11=ch0
  cp21=ch2
  IF (lb == n) cbj(lb+1)=2.0_EP*(lb+v0)/z*cbj(lb)-cbj(lb-1)
  IF (ABS(cbj(0)) > ABS(cbj(1))) THEN
     cby(lb+1)=(cbj(lb+1)*cyv0-2.0_EP*cp11/(pi*z))/cbj(0)
     cby(lb)=(cbj(lb)*cyv0+2.0_EP*cp12/(pi*z))/cbj(0)
  ELSE
     cby(lb+1)=(cbj(lb+1)*cyv1-2.0_EP*cp21/(pi*z))/cbj(1)
     cby(lb)=(cbj(lb)*cyv1+2.0_EP*cp22/(pi*z))/cbj(1)
  END IF
  cyl2=cby(lb+1)
  cyl1=cby(lb)
  DO  k=lb-1,0,-1
     ! calculate Y_{v_0+k} using backwards recurrencce for k<=B
     cylk=2.0_EP*(k+v0+1.0_EP)/z*cyl1-cyl2
     cby(k)=cylk
     cyl2=cyl1
     cyl1=cylk
  END DO
  cyl1=cby(lb)
  cyl2=cby(lb+1)
  DO  k=lb+1,n-1
     ! calculate Y_{v_0+k} using forward recurrence for k>B
     cylk=2.0_EP*(k+v0)/z*cyl2-cyl1
     cby(k+1)=cylk
     cyl1=cyl2
     cyl2=cylk
  END DO
  DO  k=2,n
     ! check if B is accurate, if not, use new B and repeat
     wa=ABS(cby(k))
     IF (wa < ABS(cby(k-1))) lb=k
  END DO
  GO TO 95
125 cdy(0)=v0/z*cby(0)-cby(1)
  DO  k=1,n
     ! Calculate the derivative of Y_v
     cdy(k)=cby(k-1)-(k+v0)/z*cby(k)
  END DO
  vm=n+v0
  RETURN
END SUBROUTINE cjyva


INTEGER FUNCTION msta1(x,mp)

  !       ===================================================
  !       Purpose: Determine the starting point for backward
  !                recurrence such that the magnitude of
  !                Jn(x) at that point is about 10^(-MP)
  !       Input :  x     --- Argument of Jn(x)
  !                MP    --- Value of magnitude
  !       Output:  MSTA1 --- Starting point
  !       ===================================================

  use constants, only : EP

  IMPLICIT real(EP) (a-h,o-z)
  real(EP), INTENT(IN OUT) :: x
  INTEGER, INTENT(IN) :: mp
  a0=ABS(x)
  n0=INT(1.1_EP*a0)+1
  f0=envj(n0,a0)-mp
  n1=n0+5
  f1=envj(n1,a0)-mp
  secant: DO  it=1,20
     ! solve eqn (5.3.4) using the secant method
     nn=int(n1-(n1-n0)/(1.0_EP-f0/f1))
     f=envj(nn,a0)-mp
     IF(ABS(nn-n1) < 1) then
        exit secant
     end IF
     n0=n1
     f0=f1
     n1=nn
     f1=f
  END DO secant
  msta1=nn
  RETURN
END FUNCTION msta1

function envj(n,x) result(evj)

  use constants, only : EP
  real(EP), intent(in) :: x
  real(EP) :: evj
  
  evj = 0.5_EP*log10(6.28_EP*n)-n*log10(1.36_EP*x/n)
  return
end function envj

INTEGER FUNCTION msta2(x,n,mp)

  !       ===================================================
  !       Purpose: Determine the starting point for backward
  !                recurrence such that all Jn(x) has MP
  !                significant digits
  !       Input :  x  --- Argument of Jn(x)
  !                n  --- Order of Jn(x)
  !                MP --- Significant digit
  !       Output:  MSTA2 --- Starting point
  !       ===================================================

  use constants, only : EP

  IMPLICIT real(EP) (a-h,o-z)
  real(EP), INTENT(IN OUT) :: x
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(IN) :: mp

  a0=ABS(x)
  hmp=0.5_EP*mp
  ejn=envj(n,a0)
  IF (ejn <= hmp) THEN
     ! if |J_N(x)|>= 10^{-p/2}, use eqn (5.3.4)
     obj=mp
     n0=INT(1.1*a0)+1
  ELSE
     ! otherwise use eqn (5.3.5)
     obj=hmp+ejn
     n0=n
  END IF
  f0=envj(n0,a0)-obj
  n1=n0+5
  f1=envj(n1,a0)-obj
  secant: DO  it=1,20
     ! solve eqn (5.3.4) or (5.3.5) using the secant method
     nn=int(n1-(n1-n0)/(1.0_EP-f0/f1))
     f=envj(nn,a0)-obj
     IF (ABS(nn-n1) < 1) then
        exit secant
     end IF
     n0=n1
     f0=f1
     n1=nn
     f1=f
  END DO secant
  msta2=nn+10
  RETURN
END FUNCTION msta2
end module complex_bessel2
