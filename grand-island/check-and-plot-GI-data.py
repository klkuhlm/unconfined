import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.mlab as mlab
from scipy.interpolate import UnivariateSpline as spline 
from glob import glob

files = sorted(glob("grand-island-test-wenzel-*.csv"))

# pumping start and end times
t0 = pylab.datestr2num('07/29/1931 06:05:00')
t1 = pylab.datestr2num('07/31/1931 06:04:00')

print 'test length in seconds:', (t1-t0)*1440*60

tlen = (t1-t0)*1440

# column 0: date 
# column 1: time
# column 2: drawdown (ft)

individualplots = True
computesplinederiv = True
drawdowncheck = True
mapcheckplot = True
saveMetricDecimalData = True

# apply uniform min/max derivatives to all data
ydmin,ydmax = (-0.2,1.0)

if drawdowncheck:
    if not individualplots:
        fig = plt.figure(2)
        ax = fig.add_subplot(111)
        colors = {'A':'red','B':'green','C':'magenta','D':'cyan',
                  'W':'pink','N':'orange','S':'black','SW':'purple'}
    
    for filename in files:
    
        if 'info' not in filename:
    
            well = filename.lstrip('grand-island-test-wenzel-').rstrip('.csv') 
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

                tv = (dt[pm]-t0)*1440
                dv = dd[pm]

                if saveMetricDecimalData:
                    # save tab-delimited file with 
                    # elapsed time in seconds
                    # and drawdown in meters

                    fh = open('out-'+filename.replace('.csv','-metric.dat'),'w')
                    for ttt,ddd in zip(tv*60,dv*0.3048):
                        fh.write('%.0f\t%.4f\n' % (ttt,ddd))
                    fh.close()

                fig = plt.figure(1)
                ax = fig.add_subplot(1,1,1)
                ax.loglog(tv,dv,'r-')
                ax.loglog(tv,dv,'k.')

                if computesplinederiv and '83' not in well:
                    ppm = np.logical_and(dt >= t0,dt <= dt[np.argmax(dd)])
                    tpv = (dt[ppm]-t0)*1440
                    dpv = dd[ppm]

                    nd = dpv.shape[0]
                    factor = nd*20
                    sp = spline(np.log(tpv),dpv,s=np.linalg.norm(dv)/factor)
                    
                    spval = np.empty((nd,2))
                    for j in range(nd):
                        spval[j,:] = sp.derivatives(np.log(tpv[j]))[0:2]
                    axd = ax.twinx()
                    ax.semilogx(tpv,spval[:,0],'b--')
                    axd.semilogx(tpv,spval[:,1],'g--')
                    axd.set_ylabel('d s/d(ln(t)) no recovery')
                    axd.set_ylim([ydmin,ydmax])

                ax.axvline(tlen)
                ax.set_ylabel('drawdown (ft)')
                if not computesplinederiv:
                    ax2 = ax.twinx()
                    ax2.semilogx(tv[1:], tv[1:]-tv[:-1],'rx')
                    ax2.set_ylabel('$\\Delta t$ (min)')
                ax.set_xlabel('time since pumping began (min)')
                ax.set_title(well)
            
                plt.savefig(filename.replace('csv','png'))
                plt.close(1)
            else:
                # most lines are one letter, except "SW"
                if '83' in well:
                    # pumping well doesn't have a "line" but is arbitrarily added 
                    # to line A for calculations
                    line = 'A'
                elif 'SW' in well:
                    line = 'SW'
                else:
                    line = well[-1:]
                ax.loglog((dt[pm]-t0)*1440,dd[pm],'-',
                          color=colors[line],linewidth=0.25)
    
    
    if not individualplots:
        ax.axvline(tlen,color='black',linewidth=0.1)
        ax.set_ylabel('drawdown (ft)')
        ax.set_xlabel('time since pumping began (min)')
        ax.set_title('July 1931 Grand Island Test USGS WSP-887 (Wenzel, 1942)')
        plt.savefig('all-GI-data.eps')
        plt.close(2)

if mapcheckplot:

    if saveMetricDecimalData:
        fh = open('grand-island-test-wenzel-metric-info.csv','w')
        fh.write('#well,z_bot(m),z_top(m),r(m),rw_obs(m)\n')
        b = 100.0*0.3048  # aquifer thickness in meters

    # angle line makes on map (east=0 & 360, north=90, west=180, south=270)
    # computed from map in Wenzel 
    angles = {'A':114,'B':205,'C':294,'D':25,
              'W':160,'SW':186,'S':238,'N':70}

    # col 0: well (number)
    # col 1: line (letter)
    # col 2: diameter (inches)
    # col 3: screen length (ft)
    # col 4: screen depth (ft BMP)
    # col 5: measuring point height (ft ALS)
    # col 6: measuring point altitude (ft AMSL)
    # col 7: distance from pumped well (ft)
    # col 8: initial water level (ft BMP)

    ddt = np.dtype([('id','S2'),('line','S2'),('screen','f8'),
                    ('bot','f8'),('measpt','f8'),
                    ('elev','f8'),('r','f8'),
                    ('wl','f8')])

    # compute x & y coordinates (pumping well is zero)
    # for creating maps for checking, etc.
    xydt = np.dtype([('id','S2'),('x','f8'),('y','f8')])

    d = np.genfromtxt('grand-island-test-wenzel-info.csv',
                      dtype=ddt,delimiter=',',
                      usecols=(0,1,3,4,5,6,7,8),skiprows=1,
                      filling_values=np.NaN)

    xy = np.empty(d.shape,dtype=xydt)    
    conv = np.pi/180.0
    xy['id'] = d['id']

    for well in d['id']:
        thiswell = xy['id'] == well
        lin = d[thiswell]['line'][0]
        #if lin == 'SW':
        #    # the ray comes from 84, not 83. Radial distance -> 83 is
        #    # correct, but XY computed is not
        #    pass
        #else:
        ang = angles[lin]*conv
        r = d[thiswell]['r'][0]
        x = r*np.cos(ang)
        y = r*np.sin(ang)
        xy['x'][thiswell] = x
        xy['y'][thiswell] = y

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

    plt.figure(4,figsize=(26,15))
    plt.subplot(122)
    plt.contourf(X,Y,ls)
    plt.plot(xy['x'],xy['y'],'k.')
    for well in d['id']:
        thiswell = d['id'] == well
        plt.annotate('%s%s' % (well,d[thiswell]['line'][0]),
                     xy=(xy[thiswell]['x'],xy[thiswell]['y']),fontsize='x-small')

        if saveMetricDecimalData:
            
            # columns
            # 0) well name
            # 1) z elevation of bottom of obs well (from bottom of aquifer - 100 ft thick)
            # 2) z elevation of top of obs well (from bottom of aquifer - 100 ft thick)
            # 3) radial distance to obs well from pumping well
            # 4) radius of observation well

            # saturated thickness assumed constant everywhere 
            # zbot = thickness - depth_of_bottom_below_water_table
            zbot = b - (d[thiswell]['bot'] - d[thiswell]['wl'])*0.3048
            ztop = zbot + d[thiswell]['screen']*0.3048 

            if d[thiswell]['screen'] < 1.0:
                rwobs = 0.25*0.3048/2  # 3" diameter piezometer
            elif d[thiswell]['screen'] > 35.0:
                rwobs = 2*0.3048/2  # 2' diameter pumping well
            else:
                rwobs = 0.3048/12/2 # 1" diameter screened well

            fh.write('%s%s,%.4f,%.4f,%.4f,%.4f\n' % 
                     (well,d[thiswell]['line'][0],zbot,ztop,
                      d[thiswell]['r']*0.3048,rwobs))

    plt.colorbar(shrink=0.5)
    plt.axis('image')
    plt.grid()

    plt.subplot(121)
    plt.contourf(X,Y,wt)
    plt.plot(xy['x'],xy['y'],'k.')
    plt.colorbar(shrink=0.5)
    plt.axis('image')
    plt.grid()
    plt.savefig('grand-island-contour-maps.eps')
    plt.close(4)

    # check screen locations and distances
    plt.figure(5,figsize=(30,8))
    #plt.plot(d['r'],d['elev']-d['bot'],'r_')
    for well in d['id']:
        thiswell = d[d['id'] == well]
        plt.annotate('%s%s' % (well,thiswell['line'][0]),
                     xy=(thiswell['r'],thiswell['elev']-thiswell['wl']),
                     fontsize='xx-small')
        if thiswell['screen'] == 0:
            # 3" piezometer (point observation)
            plt.plot(thiswell['r'],thiswell['elev']-thiswell['bot'],'k.')
        else:
            # 1" observtion well w/ 18" screened interval
            # with bottom of screen at indicated depth
            plt.plot([thiswell['r'],thiswell['r']],thiswell['elev'] -
                     np.array([thiswell['bot']-thiswell['screen'],
                               thiswell['bot']]),'k-',linewidth=1.5)
    plt.plot(d['r'],d['elev']-d['wl'],'b_')
    plt.plot(d['r'],d['elev']-d['measpt'],'k_')
    ls83 = d[d['id'] == 83]['elev']
    plt.plot([1,1],[ls83,ls83-39.5],'k-',linewidth=5)
    plt.grid()
    plt.xlim([1,d['r'].max()])
    plt.xlabel('radial distance from well 83 (ft)')
    plt.ylabel('elevation')
    plt.savefig('grand-island-screen-locations.eps')
