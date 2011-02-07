tic
   clear all;
   data = load('f505_080.dat');%Read data file into array named 'data'
   h_obs = data(:,2);%Store data in array h_obs               
   invert = true;%Set 'true' for parameter estimation (inversion)
                  %Set 'false' for forward simulation
   params0 = log([1.22536937763560e-003    528.767178292438e-003    3.76654893797454e-006 252.075063936471e-003]);%Parameters
   params=params0;
   times = data(:,1);%Measurement times
   if(invert)
       options = optimset;
       % Modify options setting
       options = optimset(options,'Display', 'off');
       options = optimset(options,'Algorithm', 'levenberg-marquardt');
       params = lsqcurvefit(@run_mod,params0,times,h_obs',[],[],options);
   end
   
   times = logspace(-3,4,100); %Time vector
   run_mod(params,times);%Run forward model

   h1=load('modified.txt'); 
   loglog(data(:,1),data(:,2),'bo','LineWidth',2);
   hold on
   loglog(h1(:,1),h1(:,2),'m-','LineWidth',2);
   h2=load('modifiedrw1.txt');
   loglog(h2(:,1),h2(:,2),'r-','LineWidth',2);
   axis([1e-3 1e4 1e-2 1]);
   xlabel('t (s)');
   ylabel('s (m)')
   legend('data','r_w = 0.0','r_w = 1.0 in.');
   toc