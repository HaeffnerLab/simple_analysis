from __future__ import division
from fitdata import *
from get_data import *
import pylab as pl

def nbar(rsb, bsb):
    r = rsb[0]/bsb[0]
    r_err = np.sqrt((bsb[1]/bsb[0])**2+(rsb[1]/rsb[0])**2) * r
    err = 1/(r-1)**2 * r_err
    return r/(1-r), err

def get_exc(fit, dataobj):
    rsb_exc = fit.param_dict[0]['Amplitude'][0]
    rsb_err = fit.param_dict[0]['Amplitude'][1]
    if rsb_exc > 1:
        print "#####################"
        print "WARNING: BAD FIT"
        print "#####################"
        return max(dataobj.data[:,1]), max(dataobj.data[:,1]*.2)
    return rsb_exc, rsb_err

def get_nbar(rsb_time, bsb_time, dataobj, do_plot=False, FitFunc=PeakFit):
    dataobj.get_data(rsb_time)
    fit_r = FitFunc(dataobj.data, show_graph=do_plot)
    rsb_exc = get_exc(fit_r)
    dataobj.get_data(bsb_time)
    fit_r = FitFunc(dataobj.data, show_graph=do_plot)
    bsb_exc = get_exc()
    

def get_nbar_list(data_dict, dataobj, rsb_index=1, do_plot=False, FitFunc=PeakFit):
    bsb_index = rsb_index + 1 % 2
    nbar_list = []
    key_list = data_dict.keys()
    try:
        key_list.sort()
    except:
       print "No integer keys for dictionary - Don't worry" 
    for key in key_list:
#        print key
        dataset = data_dict[key]
        rsb_time = dataset[rsb_index]
        dataobj.get_data(rsb_time)
        fit_r = FitFunc(dataobj.data, show_graph=do_plot)
        rsb_exc = get_exc(fit_r, dataobj)
        bsb_time = dataset[bsb_index]
        dataobj.get_data(bsb_time)
        fit = FitFunc(dataobj.data, show_graph=do_plot)
        bsb_exc = get_exc(fit, dataobj)
#        print rsb_exc
#        print bsb_exc
#        print nbar(rsb_exc,bsb_exc)
#        print "##################"
        nbar_list.append(nbar(rsb_exc,bsb_exc))
    try:
        return np.array(key_list), np.array(nbar_list)
    except:
        return (key_list), np.array(nbar_list)

def get_nbar_list_temperature(data_dict, dataobj, do_plot=False, FitFunc=PeakFit):
    nbar_list = []
    key_list = data_dict.keys()
    try:
        key_list.sort()
    except:
       print "No integer keys for dictionary - Don't worry" 
    for key in key_list:
        #print key
        dataset = data_dict[key]
        dataobj.get_data(dataset+'.dir/rsb_scan')
        fit_r = FitFunc(dataobj.data, show_graph=do_plot)
        rsb_exc = get_exc(fit_r, dataobj)
        dataobj.get_data(dataset+'.dir/bsb_scan')
        fit = FitFunc(dataobj.data, show_graph=do_plot)
        bsb_exc = get_exc(fit, dataobj)
        nbar_list.append(nbar(rsb_exc,bsb_exc))
    try:
        return np.array(key_list), np.array(nbar_list)
    except:
        return (key_list), np.array(nbar_list)


def plot_heatingrate(data_dict, filename, do_show=True):
    pl.figure(201)
    color_list = ['b','r','g','k','y','r','g','b','k','y','r',]
    fmtlist = ['s','d','o','s','d','o','s','d','o','s','d','o']
    result_dict = {}
    for key in data_dict.keys():
        x = data_dict[key][0]
        y = data_dict[key][1][:,0]
        y_err = data_dict[key][1][:,1]

        p0 = np.polyfit(x,y,1)
        fit = LinFit(np.array([x,y,y_err]).transpose(), show_graph=False)
        p1 = [0,0]
        p1[0] = fit.param_dict[0]['Slope'][0]
        p1[1] = fit.param_dict[0]['Offset'][0]
        print fit
        x0 = np.linspace(0,max(x))
        cstr = color_list.pop(0)
        fstr = fmtlist.pop(0)
        lstr = key + " heating: {0:.2f} ph/ms".format((p1[0]*1e3)) 
        pl.errorbar(x/1e3,y,y_err,fmt=fstr + cstr,label=lstr)
        pl.plot(x0/1e3,np.polyval(p0,x0),cstr)
        pl.plot(x0/1e3,np.polyval(p1,x0),cstr)
        result_dict[key] = 1e3*np.array(fit.param_dict[0]['Slope'])
    pl.xlabel('Heating time (ms)')
    pl.ylabel('nbar')
    if do_show:
        pl.legend()
        pl.show()
    if filename != None:
        pl.savefig(filename)
    return result_dict
