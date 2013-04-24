# this module implements the F.R. de Hoog, J.H. Knight, and A.N. Stokes
# numerical inverse Laplace transform algorithm.
# see "An improved method for numerical inversion of Laplace
#     transforms", SIAM J. Sci. Stat. Comp., 3, 357-366, 1982.

## an implementation of the de Hoog et al. method
## assumes proper f(p) have been computed for the p
## required for the vector of t passed to this function
## -- only one log-cycle of time should be passed at once --
## (no error checking done in this regard)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deHoog(t::Array{Float64,1}, T::Float64, fp::Array{Complex{Float64},1}, 
                alpha::Float64, tol::Float64, M::Int) 

    nt = length(t)
    np = 2*M+1
    
    # initialize Q-D table
    e = zeros(Complex{Float64},np,np)
    q = zeros(Complex{Float64},np,np-1)
    d = zeros(Complex{Float64},np)
    A = zeros(Complex{Float64},np+2,nt)
    B = ones(Complex{Float64},np+2,nt)

    z = zeros(Complex{Float64},nt)
    brem = zeros(Complex{Float64},nt)
    rem = zeros(Complex{Float64},nt)

    # Re(p) -- this is the de Hoog parameter c
    gamma = alpha - log(tol)/(2.0*T)

    q[1,1] = fp[2]/(fp[1]/2.0) # half first term
    q[2:2*M,1] = fp[3:np]./fp[2:2*M]

    # rhombus rule for filling in triangular Q-D table
    for r = 1:M
        # start with e, column 2, 1:2*M-1
        maxr = 2*(M-r)
        e[1:maxr+1,r+1] = q[2:maxr+2,r] - q[1:maxr+1,r] + e[2:maxr+2,r]
        if r != M
            # start with q, column 2, 1:2*M-2
            rq = r+1
            maxr = 2*(M-rq)+1
            q[1:maxr+1,rq] = q[2:maxr+2,rq-1] .* e[2:maxr+2,rq] ./ e[1:maxr+1,rq]
        end 
    end 
            
    # build up continued fraction coefficients
    d[1] = fp[1]/2.0 # half first term
    for r = 1:M
        d[2*r]   = -q[1,r]   # even terms
        d[2*r+1] = -e[1,r+1] # odd terms
    end

    # seed A and B for recurrence
    A[2,1:nt] = d[1]

    # base of the power series
    z[1:nt] = exp(pi*t*im/T)

    # coefficients of Pade approximation
    # using recurrence for all but last term
    for n = 3:2*M+1
        A[n,:] = A[n-1,:] + d[n-1].*A[n-2,:].*z
        B[n,:] = B[n-1,:] + d[n-1].*B[n-2,:].*z
    end

    # "improved remainder" to continued fraction
    brem[1:nt] = (1.0 + (d[2*M] - d[2*M+1]).*z)/2.0
    rem[1:nt] = -brem.*(1.0 - sqrt(1.0 + d[2*M+1].*z./brem.^2))

    # last term of recurrence using new remainder
    A[2*M+2,:] = A[2*M+1,:] + rem.*A[2*M,:]
    B[2*M+2,:] = B[2*M+1,:] + rem.*B[2*M,:]

    # diagonal Pade approximation
    # F=A/B represents accelerated trapezoid rule
    ft = exp(gamma*t)/T .* real(A[2*M+2,:]./B[2*M+2,:])
    return ft
    
end 

function pvalues(T::Float64, alpha::Float64, tol::Float64, M::Int)

    p = zeros(Complex{Float64},2*M+1)

    # real portion is constant
    sigma = alpha - log(tol)/(2.0*T)

    for i = 0:2*M
        p[i+1] = complex(sigma,pi*i/T)
    end
      
    return p

end 

# WTF. can't raise a complex number to an integer negative power?
# test function
f(p) = p^(-4.0) 

tol = 1.0E-8
alpha = 1.0E-8
M = 10

fp = zeros(Complex{Float64},2*M+1)
t = logspace(-1,1,30)
ft = zeros(Float64,length(t))

for j = 1:length(t)
    p = pvalues(t[j]*2.0,alpha,tol,M)
    for i = 1:2*M+1
        fp[i] = f(p[i])
    end
    ft[j] = deHoog(t[j:j],t[j]*2.0,fp,alpha,tol,M)[1]
    @printf("%2i %11.7f %11.7f %11.7f\n",j,t[j],t[j]^3.0/6,ft[j])
end