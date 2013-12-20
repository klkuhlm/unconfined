function spD = hl_spD(a,zD,p,params)
% Unconfined SP model based on Mishra-Neuman (2010) flow model
kappa = params(1);
beta = params(5);
u0 = beta/2.0;
d = params(8);
beta3 = params(9);

theta_hat = params(2)*exp(-params(6));
eta1 = sqrt((p*theta_hat + a.^2)/kappa);

v = sqrt(1 + (eta1/u0).^2);
v_hat = sqrt(1+(a/(u0*d)).^2);

u = u0*(1-v);
u_hat = u0*d*(v_hat-1);

chi1 = u.^2 - beta3*u - a.^2;
chi2 = (u - beta3 + beta*d).*(u - beta3) - a.^2;
chi = chi1./chi2;

if(zD>1.0)
    %Unsaturated zone
    A1 = exp(u_hat).*(phiD20(a,1.0,p,params) + (1-chi).*hl_sD(a,1.0,p,params));
    spD = A1.*exp(-u_hat*zD) + hl_phiDp(a,zD,p,params);
elseif(zD<0.0)
    %Aquitard
    spD = (phiD20(a,0.0,p,params) + hl_sD(a,0.0,p,params)).*exp(a*zD);
else
    %Saturated zone
    spD = phiD20(a,zD,p,params) + hl_sD(a,zD,p,params); 
end