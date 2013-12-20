 function spD = phiD20simp(a,zD,p,params)
%Unconfined SP model based on Mishra-Neuman (2010) flow model
kappa = params(1);%Kz/Kr, anisotropy ratio
sigD3 = params(7);
sigD1 = params(10);
beta = params(5);
theta = params(2)/beta;
aD = kappa/theta;

Delta = a.*abs(sigD3*sigD1*sinh(a) + (sigD1+sigD3)*cosh(a));

Z1 = (p/aD - a*sigD1).*hl_sDsimp(a,1.0,p,params);
Z2 = -sigD3*hl_sDsimp(a,0.0,p,params);

w1 = sigD3*sinh(a*zD) + cosh(a*zD);
w2 = sigD1*sinh(a*(1-zD)) + cosh(a*(1-zD));

spD = (Z1.*w1 + a.*Z2.*w2)./Delta;
end