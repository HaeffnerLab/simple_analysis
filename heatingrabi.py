import simple_analysis.get_data as gd
import lmfit_rabi as lmf
import numpy as np
import pylab as pl
import scipy.constants as const
import lmfit
from scipy.special.orthogonal import eval_genlaguerre as laguerre


def hrate(pars, x, data=None, eps=None):
    off = pars['offset'].value
    k = pars['rate'].value
    model = off + k * x
    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model - data)/eps



class HeatingRabi():

    def __init__(self, date, w0 = 2*np.pi*500e3, phi=np.pi/4, nmax=3000,
                 basepath='/home/viellieb/work_resonator/datavault/Experiments.dir'):
        self.date = date
        self.basepath = basepath
        k = 2 * np.pi / 729.14e-9
        eta0 = k * np.sqrt(const.hbar/(2*40*const.m_p*w0))
        self.eta = eta0 * np.cos(phi)
        self.delay = 6.7
        self.nmax = nmax
        print "Generating coupling tables"
        coupling_func = lambda n: np.exp(-1./2*self.eta**2) * laguerre(n, 0, self.eta**2)
        coupling_array = np.array([coupling_func(n) for n in range(self.nmax)])
        self.rabi_coupling = coupling_array
        self.n_array = np.arange(self.nmax)

    def Rabiflop(self, pars, t, data=None, eps=None):
        ones = np.ones_like(t)
        n = self.n_array
        omega = self.rabi_coupling
        delta = pars['delta'].value
        time_2pi = pars['time_2pi'].value
        nbar = pars['nbar'].value

        p_n = 1./ (nbar + 1.0) * (nbar / (nbar + 1.0))**n
        model = np.outer(p_n, ones) * (np.sin( np.outer( np.sqrt(omega**2+delta**2)*np.pi/time_2pi, t ))**2)
        model = np.sum(model, axis = 0)
        if data is None:
            return model
        if eps is None:
            return (model - data)
        return (model - data)/eps

    def minimize(self, data, params):
        self.result = lmfit.minimize(self.Rabiflop, params, args = (data[:,0], data[:,1], data[:,2]+.01))


    def fit_t2pi(self, fname, time_2pi_guess, nbar_guess=50, fit_nbar=True, do_plot=False):
        do = gd.ReadData(self.date,experiment='RabiFlopping', basepath=self.basepath)
        dat = do.get_data(fname)
        dat[:,0] = dat[:,0] - self.delay
        params = lmfit.Parameters()
        params.add('nbar', value=nbar_guess, vary=fit_nbar, min=0.)
        params.add('delta', value=0, min=-0.05, max=.1, vary=False)
        params.add('time_2pi', value=time_2pi_guess, vary=True, min=0.1)
        params.add('eta', value=self.eta, vary=False, min=0)
        self.minimize(dat, params)
        if do_plot:
            t = np.linspace(0,dat[:,0].max(),50)
            y = self.Rabiflop(self.result.params, t)
            pl.figure()
            pl.plot(dat[:,0],dat[:,1])
            pl.plot(t,y)
            pl.xlabel('excitation time [ms]')
            pl.ion()
        return self.result.params['time_2pi'].value

    def get_heating_rate(self, flist,tlist, time_2pi, do_plot=False, fit_t2pi=False):
        do = gd.ReadData(self.date,experiment='RabiFlopping', basepath=self.basepath)
        nlist = []
        nlist_std = []
        if do_plot:
            f1, alist = pl.subplots(len(flist))   
            alist = list(alist)
        for fname in flist:
            dat = do.get_data(fname)
            dat[:,0] = dat[:,0] - self.delay
            params = lmfit.Parameters()
            params.add('nbar', value=50, vary=True, min=0.)
            params.add('delta', value=0, min=-0.05,
                        max=.1, vary=False)
            params.add('time_2pi', value=time_2pi, vary=fit_t2pi, min=0.1)
            params.add('eta', value=self.eta, vary=False, min=0)
            t = np.linspace(0,dat[:,0].max(),100)
            self.minimize(dat, params)
#            print f.params['time_2pi']
            nlist.append(self.result.params['nbar'].value)
            nlist_std.append(self.result.params['nbar'].stderr)
            time_2pi = params['time_2pi'].value
            if do_plot:
                y = self.Rabiflop(self.result.params, t)
                a = alist.pop(0)
                a.plot(dat[:,0],dat[:,1])
                a.plot(t,y)
                a.set_yticks([0, 0.5, 1.0])
                pl.xlabel('excitation time [ms]')



        nlist = np.array(nlist)
        nlist_std = np.array(nlist_std)
        ph = lmfit.Parameters()
        ph.add('offset',value=0.1,min=0.0,vary=True)
        ph.add('rate',value=0.1,vary=True)
        #res = lmfit.minimize(hrate, ph, args=(tlist,nlist,nlist_std))
        res = lmfit.minimize(hrate, ph, args=(tlist,nlist))
        if do_plot:
            pl.figure()
            pl.errorbar(tlist,nlist,nlist_std, fmt='d')
            pl.plot(tlist,hrate(ph, tlist))
            pl.xlabel('heating time [ms]')
            pl.ylabel('nbar')
            pl.ion()
            print "Fit result (ph/ms): " + str(ph['rate'].value) + " +- " + str(ph['rate'].stderr)
            return ph['rate'].value, ph['rate'].stderr
