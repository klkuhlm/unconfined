function u = phiD(p,zD,params)
format long eng;
kappa = params(1);%Kz/Kr, anisotropy ratio
beta0 = params(10);
beta1 = params(11);
beta2 = params(12);

B1 = p*beta0*exp(-beta2)/kappa;
u = 1.0i*sqrt(4*B1/(beta1^2))*exp(beta1*zD/2);
end