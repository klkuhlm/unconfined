 function spD = phiD20(a,zD,p,params)
%Unconfined SP model based on Mishra-Neuman (2010) flow model
kappa = params(1);%Kz/Kr, anisotropy ratio
sigD3 = params(7);
beta = params(5);
d = params(8);
beta3 = params(9);

theta_hat = params(2)*exp(-params(6));

u0 = beta/2.0;
eta1 = sqrt((p*theta_hat + a.^2)/kappa);

v = sqrt(1+(eta1/u0).^2);
v_hat = sqrt(1+(a/(u0*d)).^2);

u = u0*(1-v);
u_hat = u0*d*(v_hat-1.0);

chi1 = u.^2 - beta3*u - a.^2;
chi2 = (u - beta3 + beta*d).*(u - beta3) - a.^2;
chi = chi1./chi2;

Delta = abs((u_hat*sigD3 + a).*sinh(a) + (a*sigD3 + u_hat).*cosh(a));

Z1 = ((u_hat + u).*(chi-1) - beta*d*chi).*hl_sD(a,1.0,p,params);
Z2 = sigD3*hl_sD(a,0.0,p,params);

w1 = sigD3*sinh(a*zD) + cosh(a*zD);
w2 = u_hat.*sinh(a*(zD-1)) - a.*cosh(a*(zD-1));

spD = (Z1.*w1 + Z2.*w2)./Delta;
end