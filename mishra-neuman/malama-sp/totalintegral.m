function u=totalintegral(N,M,p,rD,zD,J0,params)
%x=0.0;
x = zeros(M-N,1);
parfor k=1:M+1
    n=k-1;
    x(k) = (((-1)^n)*del_kIm(n,rD,J0,p,zD,params)/(2^(n+1)));
end
u=sum(x);