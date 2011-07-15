import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import subprocess

PIPE = subprocess.PIPE

zd = np.linspace(0.0,1.0,11) # piezometer locations
ld = np.linspace(0.005,0.995,5) # length of screen interval

c = ['red','green','blue','cyan']
l = ['dashed','dotted','dashdot']

infn =  'hantush-test-z.in'
outfn = 'hantush-test-z.out'

inputfile = """0  %i  T  T  T                :: quiet output?, model choice, dimensionless output?, timeseries?, piezometer?
2.0D-2                      :: Q (volumetric pumping rate) [L^3/T]
%.4e    %.4e               :: l/b, d/b (normalized __depth__ to bottom/top of screened interval)
2.54D-2   2.54D-2                       :: rw, rc (radius of well casing and tubing)
1.0D+0                        :: gamma (dimensionless wellbore skin)
001  0.0D+0  1.0D+0            :: pumping well time behavior and parameters 
1.00D+1                        :: b (initial saturated thickness)
1.0D-4     1.0D-1                :: Kr,kappa (radial K and ratio Kz/Kr)
1.0D-6     2.5D-1               :: Ss,Sy
1.0D+1                        ::  beta Malama linearization parameter
2.95D+0  3.7D-1   2.0D+0  2.2D-1   2.0D+1      :: Mishra/Neuman 2010; a_c,a_k,  psi_a,psi_k,  L
12  1.0D-8  1.0D-9           :: deHoog invlap;  M,alpha,tol
7  5                         :: tanh-sinh quad;  2^k-1 order, # extrapollation steps
1  1  10  50                 :: G-L quad;  min/max zero split, # zeros to integrate, # abcissa/zero
timedata.dat   15.5            :: file to read time data from (and default t otherwise)
spacedata.dat  2.5D0          :: file to read space data from (and default r otherwise)
%.4e  0.0D0  5  2.54D-2  20.0                  :: relative top obs well screen OR piezometer loc, bottom screen, quadrature order across screen (not for piezometer)
%s"""

fh = open(infn,'w')
outtheisfn = 'theis-test.out'
fh.write(inputfile % (0,1.0,0.0,0.5,outtheisfn))
fh.close()

args = ['./unconfined',infn]
print 'running theis'
stdout,stderr = subprocess.Popen(args,stdout=PIPE,stderr=PIPE).communicate()
tt,th,tdh = np.loadtxt(outtheisfn,skiprows=20,unpack=True)

fh = open(infn,'w')
outfullpenfn = 'hantush-test-fullpen.out'
fh.write(inputfile % (1,1.0,0.0,0.5,outfullpenfn))
fh.close()

args = ['./unconfined',infn]
print 'running fully-penetrating hantush'
stdout,stderr = subprocess.Popen(args,stdout=PIPE,stderr=PIPE).communicate()
ht,hh,hdh = np.loadtxt(outfullpenfn,skiprows=20,unpack=True)

for lval in ld:
    # assume screen is centered in aquifer
    ld = 0.5 + lval/2.0 # depth to bottom of screen
    dd = 0.5 - lval/2.0 # depth to top of screen

    fig = plt.figure(1,figsize=(14,12))
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
        

    ax1.loglog(tt,th,color='black',linestyle='solid',label='Theis')
    ax2.loglog(tt,tdh,color='black',linestyle='solid',label='Theis')
        

    ax1.loglog(ht,hh,color='pink',linestyle='dashed',label='Hantush full')
    ax2.loglog(ht,hdh,color='pink',linestyle='dashed',label='Hantush full')

    for i,zval in enumerate(zd):
        fh = open(infn,'w')
        fh.write(inputfile  % (1,ld,dd,zval,outfn))
        fh.close()

        args = ['./unconfined',infn]
        print 'running',ld,dd,zval
        stdout,stderr = subprocess.Popen(args,stdout=PIPE,stderr=PIPE).communicate()
        # currently not doing anything with stdout/stderr

        t,h,dh = np.loadtxt(outfn,skiprows=20,unpack=True)
        ax1.loglog(t,h, linestyle=l[i//len(c)],color=c[i%len(c)],label='$z_D=$%.1f' % zval)
        ax2.loglog(t,dh,linestyle=l[i//len(c)],color=c[i%len(c)],label='$z_D=$%.1f' % zval)
        
    ax1.set_title('partial penetration %.3f-%.3f' % (ld,dd) )
    ax2.set_xlabel('$t_D$')
    ax1.set_ylabel('$s_D$')
    ax2.set_ylabel('$d s_D/\\log(t)$')
        
    ax2.set_ylim([1.0E-5,30])
    ax1.set_ylim([1.0E-5,30])
    ax2.set_xlim([1.0E-3,1.0E2])
    ax1.set_xlim([1.0E-3,1.0E2])
        
    plt.legend(loc='lower right')
        
    plt.savefig('theis-hantush-%.3f-%.3f-compare.png' % (ld,dd))
    plt.close(1)

