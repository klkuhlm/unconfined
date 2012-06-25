import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.mlab as mlab
from glob import glob

files = sorted(glob("grand-island-test-wenzel - *.csv"))

t0 = pylab.datestr2num('07/29/1931 06:05:00')
t1 = pylab.datestr2num('07/31/1931 06:04:00')

tlen = (t1-t0)*1440

# column 0: date 
# column 1: time
# column 2: drawdown (ft)

individualplots = False
drawdowncheck = False
mapcheckplot = True

if drawdowncheck:
    if not individualplots:
        fig = plt.figure(2)
        ax = fig.add_subplot(111)
        colors = {'A':'red','B':'green','C':'magenta','D':'cyan',
                  'W':'pink','N':'orange','S':'black','SW':'purple'}
    
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
        
            if individualplots:
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
            else:
                # most lines are one letter, except "SW"
                if well[:-1].isdigit():
                    line = well[-1]
                else:
                    line = well[-2:]
                ax.loglog((dt[pm]-t0)*1440,dd[pm],'-',
                          color=colors[line],linewidth=0.25)
    
    
    if not individualplots:
        ax.set_ylabel('drawdown (ft)')
        ax.set_xlabel('time since pumping began (min)')
        plt.savefig('all-GI-data.eps')
        plt.close(2)

if mapcheckplot:

    # angle line makes on map (east=0 & 360, north=90, west=180, south=270)
    # computed from map in Wenzel 
    angles = {'A':114,'B':205,'C':294,'D':25,
              'W':160,'SW':186,'S':238,'N':70}

    # col 0: well (number)
    # col 1: line (letter)
    # col 2: diameter (inches)
    # col 3: screen depth (ft BMP)
    # col 4: measuring point height (ft ALS)
    # col 5: measuring point altitude (ft AMSL)
    # col 6: distance from pumped well (ft)
    # col 7: initial water level (ft BMP)

    fig = plt.figure(3)
    fn = glob('grand-island*info*.csv')[0]
    fh = open(fn,'rU')
    
    ddt = np.dtype([('id','S2'),('line','S2'),
                    ('srcd','f8'),('measpt','f8'),
                    ('elev','f8'),('dist','f8'),
                    ('wl','f8')])

    # compute x & y coordinates (pumping well is zero)
    # for creating maps for checking, etc.
    xydt = np.dtype([('id','S2'),('x','f8'),('y','f8')])

    d = np.genfromtxt(fh,dtype=ddt,delimiter=',',
                  usecols=(0,1,3,4,5,6,7),skiprows=1,
                  filling_values=np.NaN)

    xy = np.empty(d.shape,dtype=xydt)
    
    conv = np.pi/180.0
    xy['id'] = d['id']

    xv = []
    yv = []

    for well in d['id']:
        thiswell = xy['id'] == well
        lin = d[thiswell]['line'][0]
        ang = angles[lin]*conv
        r = d[thiswell]['dist'][0]
        x = r*np.cos(ang)
        y = r*np.sin(ang)
        xy['x'][thiswell] = x
        xy['y'][thiswell] = y

    #print xy

    buf = 50
    nx,ny = (40,40)
    x0 = xy['x'].min() - buf
    x1 = xy['x'].max() + buf
    y0 = xy['y'].min() - buf
    y1 = xy['y'].max() + buf

    X,Y = np.mgrid[x0:x1:nx*1j,y0:y1:ny*1j]

    # contour up land surface / water table
    ls = mlab.griddata(xy['x'],xy['y'],d['elev']-d['measpt'],X,Y)
    wt = mlab.griddata(xy['x'],xy['y'],d['elev']-d['wl'],X,Y)

    plt.figure(4,figsize=(17,11))
    plt.subplot(122)
    plt.contourf(X,Y,ls)
    plt.plot(xy['x'],xy['y'],'k.')
    plt.colorbar(shrink=0.5)
    plt.axis('image')
    plt.grid()

    plt.subplot(121)
    plt.contourf(X,Y,wt)
    plt.plot(xy['x'],xy['y'],'k.')
    plt.colorbar(shrink=0.5)
    plt.axis('image')
    plt.grid()
    plt.savefig('grand-island-contour-maps.png')
    plt.close(4)

