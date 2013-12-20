tic
   clear all;
 %  data = load('f505_080.dat');%Read data file into array named 'data'
  % h_obs = data(:,2);%Store data in array h_obs               
   %invert = true;%Set 'true' for parameter estimation (inversion)
                  %Set 'false' for forward simulation
   params0 = log([1e-004 1.0 1e-004 3.0e-1 5e-1 0.025 0.02 1.0 1e-3 1e-0 1e-0]);%Parameters
   params=params0;
   %times = data(:,1);%Measurement times
   %if(invert)
    %   options = optimset;
       % Modify options setting
     %  options = optimset(options,'Display', 'off');
      % options = optimset(options,'Algorithm', 'levenberg-marquardt');
       %params = lsqcurvefit(@mishraneuman,params0,times,h_obs',[],[],options);
   %end
   
   times = logspace(-2.5,6,100); %Time vector
   mishraneuman(params,times);%Run forward model

   %h0=load('mishra_ak1ac1.txt'); 
   %loglog(h0(:,1),h0(:,2),'k-','LineWidth',2);
   %hold on;
   %h0=load('mishra_ak1ac2.txt'); 
   %loglog(h0(:,1),h0(:,2),'c-','LineWidth',2);
   h0=load('mishra.txt');
   loglog(h0(:,1),h0(:,2),'m-','LineWidth',2);
   axis([1e-3 1e6 1e-1 20]);
   xlabel('t (s)');
   ylabel('s (m)')
   legend('a_{c}=1.4');
   %hold off;
   toc