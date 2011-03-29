function u=totalintegral(N,M,p,rD,zD,J0,params)
x=0.0;
for k=N:M
    x = x + (((-1)^k)*del_kIm(k,rD,J0,p,zD,params)/(2^(k+1)));
end
u=x;