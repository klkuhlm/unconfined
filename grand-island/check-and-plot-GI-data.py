import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from glob import glob

files = sorted(glob("grand-island-test-wenzel - *.csv"))

t0 = pylab.datestr2num('07/29/1931 06:05:00')
t1 = pylab.datestr2num('07/31/1931 06:04:00')

for file in files:
    print file.lstrip('grand-island-test-wenzel')
    dt,dd = np.loadtxt(file,delimiter=',',skiprows=1,
                       usecols=(2,5),unpack=True,
                       converters={2:pylab.datestr2num})

    pm = dt > t0
    rm = dt > t1

    fig = plt.figure(1)
    ax = fig.add_subplot(2,1,1)
    ax.loglog(dt[pm]-t0,dd[pm],'r-')
    ax.set_ylabel('drawdown (ft)')
    ax2 = ax.twinx()
    ax2.semilogx(dt[pm][1:]-t0, (dt[pm][1:]-dt[pm][:-1])*1440,'ro')
    ax2.set_ylabel('$\\Delta t$ (min)')
    ax.set_xlabel('time since pumping began (min)')

    ax = fig.add_subplot(2,1,2)
    ax.loglog(dt[rm]-t1,dd[rm],'g-')
    ax.set_ylabel('drawdown (ft)')
    ax2 = ax.twinx()
    ax2.semilogx(dt[rm][1:]-t1, (dt[rm][1:]-dt[rm][:-1])*1440,'ko')
    ax2.set_ylabel('$\\Delta t$ (min)')
    ax.set_xlabel('time since recovery began (min)')

    plt.savefig(file.replace('csv','png'))
    plt.close(1)
