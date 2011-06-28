import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob

fig = plt.figure()
ax = fig.add_subplot(111)

t,h,dh = np.loadtxt('theis-test.out',skiprows=20,unpack=True)
ax.loglog(t,h,label='Theis')

for fn in glob('hantush-test-*.out'):
    t,h,dh = np.loadtxt(fn,skiprows=20,unpack=True)
    ax.loglog(t,h,label=fn.replace('hantush-test-','').replace('.out',''))

ax.set_ylim([1.0E-5,30])
plt.legend()
plt.savefig('theis-hantush-compare.png')

