import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from glob import glob

files = sorted(glob("grand-island-test-wenzel - *.csv"))

t0 = pylab.datestr2num('07/29/1931 06:05:00')
t1 = pylab.datestr2num('07/31/1931 06:04:00')

tlen = (t1-t0)*1440

# column 0: date 
# column 1: time
# column 2: drawdown (ft)

for filename in files:

    if 'info' not in filename:

        well = filename.lstrip('grand-island-test-wenzel - ').rstrip('.csv') 
        print well

        t = []
        d = []

        fh = open(filename,'rU')
        lines = fh.readlines()[1:] # skip header row
        fh.close()
        
        for line in lines:
            fields = line.split(',')
            t.append(pylab.datestr2num('%s %s' % (fields[0],fields[1])))
            d.append(float(fields[2]))

        dt = np.array(t)
        dd = np.array(d)

        pm = dt > t0
        rm = dt > t1
    
        fig = plt.figure(1)
        ax = fig.add_subplot(1,1,1)
        ax.loglog((dt[pm]-t0)*1440,dd[pm],'r-')
        ax.loglog((dt[pm]-t0)*1440,dd[pm],'k.')
        ax.axvline(tlen)
        ax.set_ylabel('drawdown (ft)')
        ax2 = ax.twinx()
        ax2.loglog((dt[pm][1:]-t0)*1440, (dt[pm][1:]-dt[pm][:-1])*1440,'rx')
        ax2.set_ylabel('$\\Delta t$ (min)')
        ax.set_xlabel('time since pumping began (min)')
    
        plt.savefig(filename.replace('csv','png'))
        plt.close(1)
