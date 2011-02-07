function u=del_kIm(k,rD,J0,p,zD,params)
x=0.0;
for m = 0:k
    x = x + ((-1)^m)*nchoosek(k,m)*integral_Ik(k-m,rD,J0,p,zD,params);
end
u=x;