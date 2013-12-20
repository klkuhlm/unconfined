function u=del_kIm_sp(k,rD,J0,p,zD,params)
%x=0.0;
%%n=k
x = zeros(k,1);
parfor m = 1:k+1
    n=m-1;
    x(m) = ((-1)^n)*nchoosek(k,n)*integral_Iksp(k-n,rD,J0,p,zD,params);
end
u=sum(x);