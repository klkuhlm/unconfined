function x = hl_phiDp(a,zD,p,params)
% Particular solution for beta3 = beta*(d+1)
kappa = params(1);
beta = params(5);
d = params(8);
beta3 = params(9);

theta_hat = params(2)*exp(-params(6));
eta1 = sqrt((p*theta_hat + a.^2)/kappa);

u0 = beta/2.0;
v = sqrt(1 + (eta1/u0).^2);
u = u0*(1-v);

chi1 = u.^2 - beta3*u - a.^2;
chi2 = (u - beta3 + beta*d).*(u - beta3) - a.^2;
chi = chi1./chi2;

x = (Clr(zD,params).*chi).*hl_sD(a,zD,p,params);