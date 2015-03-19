from simple_analysis import get_data as gd
from simple_analysis import cam_postanalyze as post
import seaborn as sns
import numpy as np
import pylab as pl
from scipy import optimize

def squared_loss(params, y):
    C, phi, b = params
    x = np.array([0., 90., 180., 270.])
    e = y*(1-y)/np.sqrt(100.)
    dy = y - ((C/2)*(1 + np.sin(np.pi*x/180. + phi)) + b)
    return np.sum(0.5 * (dy / e) ** 2)

md = gd.MeasDay('2015Mar17')
rd = gd.ReadData('2015Mar17', experiment='RamseyScanGapContrast')

timestr_list = ['2028_20','2109_33','2137_27','2211_32']
c_dict = {}
for timestr in timestr_list:
    statelist = post.get_statelist(timestr, md, rd)

    statelist = np.array(statelist)
    plist = statelist.reshape([statelist.shape[0]/4,4])

    clist = []

    for i in range(plist.shape[0]):
        p = plist[i,:]
        c, phase_shift, b = optimize.fmin(squared_loss, [0,0,0], args=tuple([p]) )
        clist.append(np.abs(c))
    c_dict[timestr] = clist
    pl.figure()
    for i in range(4):
        pl.plot(plist[:,i])

pl.figure()
for i, timestr in enumerate(timestr_list):
    clist = c_dict[timestr]
    pl.subplot(4,1,i)
    pl.plot(clist, label=timestr)
    pl.legend()
pl.show()
