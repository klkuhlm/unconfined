function u = Clr(zD,params)
d = params(8);
beta = 0*params(5);
%phi_Dak = params(6);
beta3 = 0.0*params(9);

if(zD>1)
    %u = exp(phi_Dak - beta*d*(zD-1));
    u = exp(-(beta3 - beta*d)*(zD-1));
else
    u = 1.0;
end