tic
   clear all;
 %  data = load('f505_080.dat');%Read data file into array named 'data'
  % h_obs = data(:,2);%Store data in array h_obs               
   %invert = true;%Set 'true' for parameter estimation (inversion)
                  %Set 'false' for forward simulation
   params0 = log([1.22536937763560e-003    0.528767178292438    3.76654893797454e-006 0.2521]);%Parameters
   params=params0;
   %times = data(:,1);%Measurement times
   %if(invert)
    %   options = optimset;
       % Modify options setting
     %  options = optimset(options,'Display', 'off');
      % options = optimset(options,'Algorithm', 'levenberg-marquardt');
       %params = lsqcurvefit(@run_mod,params0,times,h_obs',[],[],options);
   %end
   
   times = logspace(-3,4,100); %Time vector
   run_mod(params,times);%Run forward model

%    h1=load('mishra_ac1ak1.txt'); 
%    loglog(h1(:,1),h1(:,2),'b-','LineWidth',2);
%    hold on;
%    h2=load('mishra_ac1ak2.txt'); 
%    loglog(h2(:,1),h2(:,2),'r-','LineWidth',2);
%    h3=load('mishra_ac1ak3.txt'); 
%    loglog(h3(:,1),h3(:,2),'m-','LineWidth',2);
%    h4=load('mishra.txt'); 
%    loglog(h4(:,1),h4(:,2),'m-','LineWidth',2);
%    axis([1e-3 1e4 1e-2 0.5]);
%    xlabel('t (s)');
%    ylabel('s (m)')
%    legend('a_k = 20.4','a_k = 15.4','a_k = 10.4','a_k = 5.4');
%    hold off;
   toc