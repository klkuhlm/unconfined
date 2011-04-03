function s = run_mod(params0,times)
fid = fopen('mishra.txt', 'w');
%format long eng;
Nk = 8; %Number of terms in Laplace inversion Fourier series
N=0; %First zero of BesselJ0 at which to start acceleration
M=10; %Number of zeros of BesselJ0 to use in Hankel transform inversion
J0=load('besJ0zeros.dat');
estimates = exp(params0)

Kr = estimates(1); %Horiz hydraulic conductivity
kappa = estimates(2); %anisotropy
Ss = estimates(3); %Specific storage
Sy = estimates(4); %Specific yield

Q = 320.0; %pumping rate in gpm
rw = 0.0254; %Pumping well radius
rw_obs = 0.0254; %Observation well radius
b = 52.669; %Saturated thickness
d = 0.0234; %Depth to top of test interval
l = 50.288; %Depth to bottom of test interval
z = 34.564; %Vertical position of observation point
r = 6.5837; %Radial position of observation point
Lc = b;
sF = 20.0; %Observation well shape factor
L = 20.0; %Thickness of unsaturated zone
ak = 5.4; %Unsaturated conductivity exponential decay constant
ac = 25.2; %Effective saturation exponential decay constant
phi_a = 0.739; %Air entry pressure head
phi_k = 0.521; %Pressure head above which k0 = 1.0

bD = b/Lc;
dD = d/Lc;
lD = l/Lc;
rDw = rw/Lc;
rDw_obs = rw_obs/Lc;
rD = r/Lc;
zD = z/Lc;

alpha = Kr/Ss;
sig = Sy/(Lc*Ss);
Tc = (Lc^2)/alpha;
Hc = Q*6.309e-5/(4*pi*b*Kr);

params(1) = kappa;
params(2) = sig;
params(3) = bD;
params(4) = dD;
params(5) = lD;
params(6) = rDw;
params(7)= rDw^2/(2*(l-d)*Ss);
params(8) = 0.3; %\beta_D, the linearization parameter of Malama (2011)
                 %\beta_D = 0 corresponds to Neuman (1972, 1974)
params(9) = pi*rDw_obs^2/(sF*Ss);%Dimensionless delay time for obs well
params(10) = ac*Sy/Ss;
params(11) = b*(ak-ac);
params(12) = ak*(phi_a - phi_k);
params(13) = ak*b;
params(14) = L/b;

tD = 60*times/Tc;
Ntimes = length(tD);
sD = zeros(1,Ntimes);
for n=1:Ntimes
    sD(n) = dehoog2(tD(n),rD,zD,1.2*tD(n),Nk,N,M,J0,params);
    %sD(n) = lap_invert(N,M,tD(n),rD,zD,J0, params);
    fprintf(fid, '%e %e\n',tD(n)*Tc/60,Hc*sD(n));
end
s = Hc*sD;
fclose(fid);