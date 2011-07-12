import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

t,h,dh = np.loadtxt('theis-test.out',skiprows=20,unpack=True)
ax.loglog(t,h,color='black',linestyle='solid',label='Theis')

t,h,dh = np.loadtxt('hantush-test-fullpen.out',skiprows=20,unpack=True)
ax.loglog(t,h,color='black',linestyle='dashed',label='Hantush full')

c = ['red','green','blue','cyan']
l = ['dashed','dotted','dashdot']

for i in range(11):
    if i < 10:
        zd = i/10.0 + 0.025
    else:
        zd = i/10.0 - 0.025

    fn = 'hantush-test-%.3f.out' % zd
    t,h,dh = np.loadtxt(fn,skiprows=20,unpack=True)
    ax.loglog(t,h,linestyle=l[i//len(c)],color=c[i%len(c)],label='%.3f' % zd)

ax.set_ylim([1.0E-5,30])
plt.legend(loc='lower right')
plt.savefig('theis-hantush-compare.png')

