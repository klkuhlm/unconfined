function u = hl_sDH(a,p,zD,params)
%Solution for partially penetrating pumping well
%Gives point drawdown (piezometer)
%Based on Malama et al. (2010)
%Equivalent to Hantush (1964)
%Implementation here accounts for
% (a) Pumping wellbore storage after Papadopulos & Cooper (1968)
% (b) Observation wellbore storage based on Black & Kipp (1981)
format long eng;
kappa = params(1);%Kz/Kr, anisotropy ratio
dD = params(4); %Dimensionless depth to top of pumping well
lD = params(5); %Dimensionless depth to bottom of pumping well
rDw = params(6); %Dimensionless pumping well radius
CDw = params(7); %Dimensionless pumping well storage coefficient
tDb = params(9); %Dimensionless observation well delay time

xDw = rDw*sqrt(p);
A0 = 2/(p*CDw*besselk(0,xDw) + xDw*besselk(1,xDw));
eta = sqrt((p+a.^2)/kappa);

uDf = A0./(p*(p + a.^2)*(p*tDb + 1.0));
u = (uDf/(lD-dD)).*uDp(eta,zD,params);