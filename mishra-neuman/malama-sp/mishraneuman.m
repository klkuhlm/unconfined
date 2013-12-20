function s = mishraneuman(params0,times)
%Implemetation of Mishra-Neuman (2010) unconfined flow model
%Fully penetrating pumping well, no wellbore storage
%Gives point drawdown (piezometer) in saturated & unsaturated zones
%Semi-infinte unsaturated zone and confining unit

savefile='mishra.txt';
Nk = 12; %Number of terms in Laplace inversion Fourier series
N=0; %First zero of BesselJ0 at which to start acceleration
M=6; %Number of zeros of BesselJ0 to use in Hankel transform inversion
J0=load('besJ0zeros.dat');
estimates = exp(params0)

Kr = estimates(1); %Vertical saturated zone hydraulic conductivity
kappa = estimates(2); %Horizontal hydraulic conductivity
Ss = estimates(3); %Specific storage
Sy = estimates(4); %Specific yield

Q = 10.0; %pumping rate in gpm
b1 = 10.0; %Unsaturated zone thickness
b2 = 20.0; %Saturated zone thickness
b3 = 10.0; %Confining unit thickness
z = 10.0; %Vertical position of observation point
%z = linspace(0,40,100);
r = 5.0; %Radial position of observation point

rD = r/b2;%Dimensionless radial dtistance to observation point
zD = z/b2;%Dimensionless vertical position of observation point

ac = estimates(5);%Moisture retention exponent
phi_a = estimates(6);%Air entry pressure head
phi_k = estimates(7);%Pressure head above which k0 = 1.0

alpha = Kr/Ss;%Hydraulic diffusivity
Tc = (b2^2)/alpha;%Characteristic time
Hc = Q*6.309e-5/(4*pi*b2*Kr);%Characteristic head

params(1) = kappa;%Anisotropy ratio (=K_z/K_r)
params(2) = ac*Sy/Ss;%theta, dimensionless storage ratio
params(3) = b1/b2;%bD1, dimensionless unsaturated zone thickness
params(4) = b3/b2;%bD3, dimensionless confining unit thickness
params(5) = ac*b2;%beta0, exponent in moisture retention curve
params(6) = ac*(phi_a - phi_k);%Dimensionless exponent in unsaturated hydraulic conductivity

%%matlabpool open local 8;
tD = 60*times/Tc;%Dimensionless time; note input vector times[] is in minutes
Ntimes = length(tD);%Length of input times[] vector
sD = zeros(1,Ntimes);%Initiating sD vector with zeros
for n=1:Ntimes
    %Computation of phiD using de Hoog (1980) algorithm
    sD(n) = dehoog2b(tD(n),rD,zD,1.2*tD(n),Nk,N,M,J0,params);
end
%parfor n=1:length(zD)
    %Computation of phiD using de Hoog (1980) algorithm
 %   sD(n) = dehoog2b(tD,rD,zD(n),1.2*tD,Nk,N,M,J0,params);
%end

data = [tD;sD]';%Dimensionless time & dimensionless drawdown output
save(savefile, 'data','-ascii', '-tabs');%Saving result to output file
s = Hc*sD;
%matlabpool close;
%fclose(fid);