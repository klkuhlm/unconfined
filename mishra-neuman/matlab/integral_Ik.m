function u=integral_Ik(k,rD,J0,p,zD,params)

format long;
a = J0(k+1)/rD;
b = J0(k+2)/rD;
F =@(x)(x.*(hl_sDH(x,p,zD,params)+hl_sDU(x,p,zD,params)).*abs(besselj(0,x*rD)));
u = quadgk(F,a,b);