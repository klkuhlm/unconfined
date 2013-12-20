function sD = hl_sDsimp(a,zD,p,params)
%Mishra & Neuman 2010 Unconfined flow solution
%Fully penetrating pumping well, no wellbore storage
%Gives point drawdown (piezometer)
kappa = params(1);%Kz/Kr, anisotropy ratio
beta = params(5);
theta = params(2)/beta;

eta = sqrt((p+a.^2)/kappa);
aD = kappa/theta;
Delta = (eta*aD/p).*sinh(eta) + cosh(eta);

%Saturated zone flow solution
sD = (2./(kappa*p.*eta.^2)).*(1.0 - cosh(eta*zD)./Delta);