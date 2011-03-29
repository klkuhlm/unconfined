function u = hl_sDU(a,p,zD,params)
%Solution for partially penetrating pumping well
%Gives point drawdown (piezometer)
%Based on Malama et al. (2010)
%Equivalent to Hantush (1964)
%Implementation here accounts for
% (a) Pumping wellbore storage after Papadopulos & Cooper (1968)
% (b) Observation wellbore storage based on Black & Kipp (1981)
format long eng;
kappa = params(1);%Kz/Kr, anisotropy ratio
beta1 = params(11);
beta3 = params(13);
LD = params(14);

eta = sqrt((p+a.^2)/kappa);
B2 = (a.^2)/kappa;
nu = sqrt((beta3^2 + 4*B2)/beta1^2);

v1 = phiD(p,LD,params);
v2 = phiD(p,0,params);
v3 = (beta3 + nu*beta1);

a11 = v3.*besselj(nu,v1,1) - v1*beta1*besselj(nu+1,v1,1);
a12 = v3.*bessely(nu,v1,1) - v1*beta1*bessely(nu+1,v1,1);
a21 = v3.*besselj(nu,v2) - v2*beta1*besselj(nu+1,v2);
a22 = v3.*bessely(nu,v2) - v2*beta1*bessely(nu+1,v2);

A = a12./a11;

Delta2 = a22 - A.*a21;
Delta1 = (bessely(nu,v2) - A.*besselj(nu,v2)).*(2*eta.*sinh(eta)./Delta2) - cosh(eta);

u = (hl_sDH(a,p,1.0,params)./Delta1).*cosh(eta*zD);