tic
clear all;
% params (used in later routines)
% 1   kappa (Kz/Kr)
% 2   theta (Sy*ac/Ss)
% 3   bD1 (b1/b2) dimensionless unsat zone thickness
% 4   bD3 (b3/b2) dimensionless confining unit thickness
% 5   beta0 (ac*b2) moisture retention curve exponent
% 6   ac*(phia - phik) hydraulic conductivity curve exponent
% 7
% 8
% 9
% 10 
% 11

% params used here
%              Kr      kappa    Ss      Sy      ac     phia     phik
params = log([1.0E-4, 1.0E+0, 1.0E-4, 3.0E-1, 5.0E-1, 2.5E-2, 2.0E-2, ...
	      1.0E+0, 1.0E-3, 1.0E+0, 1.0E+0]);
%              8         9       10    11  
   

times = logspace(-2.5,6,100); %Time vector
mishraneuman(params,times);%Run forward model

h0=load('mishra.txt');
loglog(h0(:,1),h0(:,2),'m-','LineWidth',2);
axis([1e-3 1e6 1e-1 20]);
xlabel('t (s)');
ylabel('s (m)')
legend('a_{c}=1.4');
toc