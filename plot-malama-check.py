import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess

PIPE = subprocess.PIPE

zdvec = np.linspace(0.0,1.0,5) # piezometer locations
sigmavec = np.logspace(np.log10(2.5E+1),np.log10(2.5E+5),5)  # Sy/(Ss*b)
kappavec = np.logspace(-3,1,5)  # Kz/Kr

c = ['red','green','blue','cyan','magenta','orange']
l = ['dashed','dotted','dashdot']

infn =  'malama-test-z.in'

inputfile = """0  %i  T  T  T                :: quiet output?, model choice, dimensionless output?, timeseries?, piezometer?
2.0D-2                      :: Q (volumetric pumping rate) [L^3/T]
6.0D-1    4.0D-1               :: l/b, d/b (normalized __depth__ to bottom/top of screened interval)
2.54D-2   2.54D-2              :: rw, rc (radius of well casing and tubing)
1.0D+0                        :: gamma (dimensionless wellbore skin)
001  0.0D+0  1.0D+0            :: pumping well time behavior and parameters 
1.00D+1                        :: b (initial saturated thickness)
1.0D-4     %.4e                :: Kr,kappa (radial K and ratio Kz/Kr)
%.4e       %.4e               :: Ss,Sy
%.4e                        ::  beta Malama linearization parameter
2.95D+0  3.7D-1   2.0D+0  2.2D-1   2.0D+1      :: Mishra/Neuman 2010; a_c,a_k,  psi_a,psi_k,  L
10  1.0D-8  1.0D-9           :: deHoog invlap;  M,alpha,tol
7  5                         :: tanh-sinh quad;  2^k-1 order, # extrapollation steps
1  1  10  50                 :: G-L quad;  min/max zero split, # zeros to integrate, # abcissa/zero
timedata.dat   15.5            :: file to read time data from (and default t otherwise)
spacedata.dat  2.5D0          :: file to read space data from (and default r otherwise)
%.4e  0.0D0  5  2.54D-2  20.0                  :: relative top obs well screen OR piezometer loc, bottom screen, quadrature order across screen (not for piezometer)
%s"""

# **************************************************

bval = 10.0  # assume thickness of 10.0
Sy = 0.25  # assume specific yeild of 25%
beta = 0.0

# **************************************************

for sigma in sigmavec:
    for kappa in kappavec:

        fig = plt.figure(1,figsize=(14,12))
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)

        for i,zval in enumerate(zdvec):

            fh = open(infn,'w')
            outhantfn = 'early-hantush-test.out'
            mod = 1
            Ss = Sy/(sigma*bval)
            fh.write(inputfile % (mod,kappa,Ss,Sy,beta,zval,outhantfn))
            fh.close()
        
            args = ['./unconfined',infn]
            print 'running early hantush',zval
            stdout,stderr = subprocess.Popen(args,stdout=PIPE,stderr=PIPE).communicate()
            t,h,dh = np.loadtxt(outhantfn,skiprows=20,unpack=True)

            ax1.loglog(t,h, linestyle='solid',color=c[i%len(c)],label='Hantush $z_D=$ %.1f' %zval)
            ax2.loglog(t,dh,linestyle='solid',color=c[i%len(c)],label='Hantush $z_D=$ %.1f' %zval)

##%##            # **************************************************
##%##            fh = open(infn,'w')
##%##            outhantfn = 'late-hantush-test.out'
##%##            mod = 1
##%##            beta = 0.0
##%##            fh.write(inputfile % (mod,kappa,Sy,Sy,beta,zval,outhantfn))
##%##            fh.close()
##%##
##%##            print 'running late hantush',zval
##%##            stdout,stderr = subprocess.Popen(args,stdout=PIPE,stderr=PIPE).communicate()
##%##            t,lh,ldh = np.loadtxt(outhantfn,skiprows=20,unpack=True)
##%##        
##%##            ax1.loglog(t,h, linestyle='dashdot',color=c[i%len(c)],label='late H $z_D=$ %.1f' %zval)
##%##            ax2.loglog(t,dh,linestyle='dashdot',color=c[i%len(c)],label='late H $s_D=$ %.1f' %zval)

            fh = open(infn,'w')
            outneumanfn = 'neuman-test.out'
            mod = 5
            beta = 0.0
            fh.write(inputfile % (mod,kappa,Ss,Sy,beta,zval,outneumanfn))
            fh.close()

            print 'running',sigma,kappa,zval
            stdout,stderr = subprocess.Popen(args,stdout=PIPE,stderr=PIPE).communicate()

            t,h,dh = np.loadtxt(outneumanfn,skiprows=20,unpack=True)
            ax1.loglog(t,h, linestyle=l[i//len(c)],color=c[i%len(c)],label='$z_D=$%.1f' % zval)
            ax2.loglog(t,dh,linestyle=l[i//len(c)],color=c[i%len(c)],label='$z_D=$%.1f' % zval)
    
        ax1.set_title('$\\sigma=$%.2e $\\kappa=$%.2e' % (sigma,kappa))
        ax2.set_xlabel('$t_D$')
        ax1.set_ylabel('$s_D$')
        ax2.set_ylabel('$d s_D/\\log(t)$')
    
        ax2.set_ylim([3.0E-6,300])
        ax1.set_ylim([3.0E-6,300])
#        ax2.set_xlim([1.0E-4,1.0E+8])
#        ax1.set_xlim([1.0E-4,1.0E+8])
    
        plt.legend(loc='lower right')
    
        plt.savefig('neuman-%.2e-%.2e-compare.png' % (sigma,kappa))
        plt.close(1)

