function sD = hl_sDavg(a,zD,p,params)
%Mishra & Neuman 2010 Unconfined flow solution
%Fully penetrating pumping well, no wellbore storage
%Gives point drawdown (piezometer)
kappa = params(1);%Kz/Kr, anisotropy ratio
beta = params(5);
bD1 = params(3);

theta_hat = params(2)*exp(-params(6));

eta = sqrt((p+a.^2)/kappa);
eta1 = sqrt((p*theta_hat + a.^2)/kappa);

u0 = beta/2.0;
v = sqrt(1+(eta1/u0).^2);
u = u0*(1-v);

Delta = (eta.*sinh(eta) - u.*cosh(eta));

if(zD<1.0)   
    %Saturated zone flow solution
    sD = (2./(kappa*p.*eta.^2)).*(1.0 + u.*sinh(eta)./(eta.*Delta));
else
    %Unsaturated zone flow solution
    A12 = 2*sinh(eta)./(p*kappa*eta.*Delta);
    sD = (A12./(u*bD1)).*(exp(u*bD1) - 1.0);
end