function u=integral_Iksp(k,rD,J0,p,zD,params)
%for solution without unsaturated flow use hl_spDsimp(x,zD,p,params)
%for solution with unsaturated flow use hl_spD(x,zD,p,params)
%format long;
a = J0(k+1)/rD;
b = J0(k+2)/rD;
F =@(x)(x.*hl_spDsimp(x,zD,p,params).*abs(besselj(0,x*rD)));
u = quadgk(F,a,b);