function Hnorm = dehoog2sp(tD,rD,zD,tDmax,Nk,N,M,J0,params)

%dehoog2: Function for calculating the inverse laplace transform needed to
%evaluate slug tests using an efficient numerical method. Original source
%is de Hoog et al. (1982) "An improved method for numerical inversion of
%Laplace transforms." This function in particular is used to calculate the
%numerical inverse laplace transform of the modified solution to Hyder et
%al. (1994) given by Malama et al. (2010) "Modeling slug tests in
%unconfined aquifers taking into account source well skin and inertial
%effects". The function calls the sub-function hyder, which is the solution
%in the laplace domain
%
%Usage:
%[Hnorm] = dehoog2(t,tmax,params,testconstants)
%where:
%   -Hnorm (scalar) is the calculated normalized head (H/H_0) for the
%   specified time
%   -t (scalar) is the time at which to compute the solution
%   -tmax (scalar) is a "late-time" defining parameter used to estimate the
%   error in the solution
%   -params (6 x 1) is a vector of aquifer parameters to be tried for
%     test-fitting. params should be formatted as follows:
%     params = [ln(Kr2) ln(Kz2) ln(Ss2) ln(Kr1) ln(Kz1) ln(Ss1)] where zone
%     2 is the main aquifer and zone 1 is the wellbore skin.
%   -testconstants (6 x 1) is a vector of constants that apply to the given
%     test. testconstants should be formatted as follows:
%     testconstants = [l d r_w r_sk r_c B]
%
%Solution by de Hoog et al., Hyder et al. and Malama et al.
%Original coding by B. Malama
%Coding / optimization by M. Cardiff

%Nk = 20;
%tmax = 20.0;
Err = 1e-6;
T = 0.8*tDmax;
beta = -log(Err)/(2*T);
alpha = 1i*pi*tD/T;
z = exp(alpha);
lapf = zeros(Nk+1,1);
for k=1:Nk+1
    p = beta + 1i*(k-1)*pi/T;
    lapf(k) = totalintegral_sp(N,M,p,rD,zD,J0,params);
end
M = Nk/2;
eps1 = zeros(Nk+1,1) + 1i*zeros(Nk+1,1);
a = lapf; a(1) = lapf(1)/2.0;

q1 = zeros(Nk,1);
q1(1:1:Nk) = a(2:1:(Nk+1))./a(1:1:Nk);

d(1) = a(1);
d(2) = -q1(1);
eps2 = zeros(Nk-3,1);
q2 = zeros(Nk-4,1);
for r=2:M+1
    for n=1:(Nk-2*r+1)
        eps2(n)=q1(n+1)-q1(n)+eps1(n+1);
    end
    for n=1:(Nk-2*r)
        q2(n)=q1(n+1)*eps2(n+1)/eps2(n);
    end
    d(2*r-1) = -eps2(1);
    if(r<M+1)
        d(2*r) = -q2(1);
    end
    eps1(1:1:(Nk-2*r+1)) = eps2(1:1:(Nk-2*r+1));
    q1(1:1:(Nk-2*r)) = q2(1:1:(Nk-2*r));
end
A=zeros(Nk+2,1) + 1i*zeros(Nk+2,1);
A(2)=d(1);
B=ones(Nk+2,1) + 1i*zeros(Nk+2,1);
for n=3:(Nk+2)
    A(n)=A(n-1)+d(n-1)*z*A(n-2);
    B(n)=B(n-1)+d(n-1)*z*B(n-2);
end

Hnorm = exp(beta*tD)*(real(A(Nk+2)/B(Nk+2)))/T;