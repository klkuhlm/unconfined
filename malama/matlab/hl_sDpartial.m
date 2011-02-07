function u = hl_sDpartial(a,p,zD,params)
format long eng;
kappa = params(1);
sig = params(2);
dD = params(4);
lD = params(5);
rDw = params(6);
CDw = params(7);
betaD = params(8);

xDw = rDw*sqrt(p);
A0 = 2/(p*CDw*besselk(0,xDw) + xDw*besselk(1,xDw));
alphaD = kappa/sig;
eta = sqrt((p+a.^2)/kappa);
xi = eta*alphaD/p;
Del = xi.*sinh(eta) + (1+betaD*eta.*xi).*cosh(eta);

uDf = A0./(p*(p + a.^2));
u = (uDf/(lD-dD)).*(uDp(a,p,zD,params) - uDp(a,p,1.0,params).*cosh(eta*zD)./Del);