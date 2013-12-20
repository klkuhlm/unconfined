function spD = hl_spDsimp(a,zD,p,params)
% Unconfined SP model based on Mishra-Neuman (2010) flow model

if(zD>1.0)
    %Unsaturated zone
    A1 = phiD20simp(a,1.0,p,params) + hl_sDsimp(a,1.0,p,params);
    spD = A1.*exp(-a*(zD-1.0));
elseif(zD<0.0)
    %Aquitard
    spD = (phiD20simp(a,0.0,p,params) + hl_sDsimp(a,0.0,p,params)).*exp(a*zD);
else
    %Saturated zone
    spD = phiD20simp(a,zD,p,params) + hl_sDsimp(a,zD,p,params); 
end