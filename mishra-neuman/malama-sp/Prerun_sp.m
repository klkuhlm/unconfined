tic
   clear all;
 %  data = load('f505_080.dat');%Read data file into array named 'data'
  % h_obs = data(:,2);%Store data in array h_obs               
   %invert = true;%Set 'true' for parameter estimation (inversion)
                  %Set 'false' for forward simulation
   params0 = log([1e-004 1.0 1e-004 3.0e-1 5e-1 0.025 0.02 1.0 1e-3 1e-0 1e-0]); %Parameters
   params=params0;
   %times = data(:,1);%Measurement times
   %if(invert)
    %   options = optimset;
       % Modify options setting
     %  options = optimset(options,'Display', 'off');
      % options = optimset(options,'Algorithm', 'levenberg-marquardt');
       %params = lsqcurvefit(@mishraneuman,params0,times,h_obs',[],[],options);
   %end
   
   times = logspace(-3,8,200); %Time vector
   %times=1e4;
   %mishraneuman(params,times);%Run forward model
   
   %h0=load('mishra.txt');
   %plot(h0(:,1),h0(:,2),'m--','LineWidth',2);
   %loglog(h0(:,1),h0(:,2),'m--','LineWidth',2);
   %axis([1e-3 1e10 1e-4 2e1]);
   %hold on;
   
   unsatsp(params,times);%Run forward model
   sp0=load('sp.txt');
   semilogx(sp0(:,1),sp0(:,2),'m-','LineWidth',2);
   hold on;
  
   %axis([1e-3 1e8 1e-4 10]);
   xlabel('t (s)');
   ylabel('\phi_{i} (V)')
   %legend('\sigma_{D,3} = 1e4');
   %hold off;
   toc