tic
   clear all;
   z = linspace(20,30,20);
   k_r = k_rel(z);
   plot(z,k_r,'r-','LineWidth',2);
   xlabel('z');
   ylabel('k_r')
   %legend('a_c=1.0','z=a_c=2.0');
   %hold off;
   toc