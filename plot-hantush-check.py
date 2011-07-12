import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(14,12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

t,h,dh = np.loadtxt('theis-test.out',skiprows=20,unpack=True)
ax1.loglog(t,h,color='black',linestyle='solid',label='Theis')
ax2.loglog(t,dh,color='black',linestyle='solid',label='Theis')

t,h,dh = np.loadtxt('hantush-test-fullpen.out',skiprows=20,unpack=True)
ax1.loglog(t,h,color='pink',linestyle='dashed',label='Hantush full')
ax2.loglog(t,dh,color='pink',linestyle='dashed',label='Hantush full')

c = ['red','green','blue','cyan']
l = ['dashed','dotted','dashdot']

for i in range(11):
    zd = i/10.0

    fn = 'hantush-test-%.1f.out' % zd
    t,h,dh = np.loadtxt(fn,skiprows=20,unpack=True)
    ax1.loglog(t,h,linestyle=l[i//len(c)],color=c[i%len(c)],label='%.1f' % zd)
    ax2.loglog(t,dh,linestyle=l[i//len(c)],color=c[i%len(c)],label='%.1f' % zd)

ax1.set_title('partial penetration 0.45-0.55')
ax2.set_xlabel('$t_D$')
ax1.set_ylabel('$s_D$')
ax2.set_ylabel('$d s_D/\\log(t)$')

ax2.set_ylim([1.0E-5,30])
ax1.set_ylim([1.0E-5,30])
ax2.set_xlim([1.0E-3,1.0E2])
ax1.set_xlim([1.0E-3,1.0E2])

plt.legend(loc='lower right')

plt.savefig('theis-hantush-0.45-0.55-compare.png')

