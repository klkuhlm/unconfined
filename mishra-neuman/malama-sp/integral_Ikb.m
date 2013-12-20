function u=integral_Ikb(k,rD,J0,p,zD,params)
%use hl_sDavg(x,zD,p,params) for vertically averaged drawdown
%use hl_sD(x,zD,p,params) for point drawdown
%format long;
a = J0(k+1)/rD;
b = J0(k+2)/rD;
F =@(x)(x.*hl_sD(x,zD,p,params).*abs(besselj(0,x*rD)));
u = quadgk(F,a,b);