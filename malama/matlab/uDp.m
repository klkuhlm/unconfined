function uD = uDp(a,p,zD,params)
format long eng;
kappa = params(1);
dD = params(4);
lD = params(5);
eta = sqrt((p+a.^2)/kappa);

f1 = sinh(eta*dD);
f2 = sinh(eta*(1-lD));
f3 = exp(-eta*(1-lD)) - (f1 + exp(-eta).*f2)./sinh(eta);
g1 = cosh(eta*(1-dD-zD));
g2 = (f1.*cosh(eta*zD) + f2.*cosh(eta*(1-zD)))./sinh(eta);
if(zD>1-dD)
    uD = g1 - g2;
elseif(zD<(1-lD))
    uD = f3.*cosh(eta*zD);
else
    uD = 1 - g2;
end