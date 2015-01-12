#!/usr/bin/env python
# -*- mode: Python; coding: latin-1 -*-
# Time-stamp: "2011-01-11 19:11:46 c704252"

#  file       fit_parity_bayesian.py
#  copyright  (c) Thomas Monz 2010

import os

import numpy as np
import numpy.matlib as npml
import scipy.special as spec
# 2011-01-11, TM: cluster only need to calculate, not plot ...
try:
    import pylab as pl
except:
    pass


def fit_contrast_bayesian(phase,parity, NumberOfIons = 2, NumberOfCycles = 100, mean_amp = 0.9, pm_mean_amp = 0.10):
    # the code is already tons better and only makes a bayesian around the value that we think it is like -> amp_mean
    gridsize_z = [200]
    gridsize_phi = [100]

    z = np.linspace(np.clip(mean_amp-pm_mean_amp,0.001,1-0.001),np.clip(mean_amp+pm_mean_amp,0.001,1-0.001),gridsize_z[0])
    phi = np.linspace(0,2*np.pi, gridsize_phi[0])
    output = np.zeros((gridsize_z[0],gridsize_phi[0]))

    # for easy of fitting we'll adjust parameters/units
    phase = np.pi*phase*NumberOfIons

    for k in xrange(len(z)):
        for l in xrange(len(phi)):
            output[k,l] = lnprobab(z[k],phi[l],phase, parity*NumberOfCycles, NumberOfCycles)

    # first lets get real probabilities
    # we don't care about the phase, so we trace over that
    # renormalise
    # find mean and std

    lnprob = output
    # multiply so that it calculatable. max value will be set to be at something like exp(50)
    lnprob -= lnprob.max() + 50

    tmp = np.exp(lnprob)
    prob = 1.*tmp/np.sum(tmp)


    # sum over phi to get z_probs
    z_probs  = np.sum(prob, axis = 1)
    z_mean = np.sum(z_probs*z)
    sigma_z = np.sqrt(np.sum(z_probs*(z-z_mean)**2))

    # sum over phi to get z_probs
    phi_probs  = np.sum(prob, axis = 0)
    phi_mean = np.sum(phi_probs*phi)
    sigma_phi = np.sqrt(np.sum(phi_probs*(phi-phi_mean)**2))
#    pl.imshow(prob)
#    pl.show()
    return (z_mean, sigma_z, phi_mean, z_probs, z)


def lnprobab(z, phi, phase, counts, NumberOfCycles):
    probk = np.zeros(len(counts))
    for k in xrange(len(counts)):
        if z <=1 and z>=0:
            probk[k] = lnbeta((1+z*np.cos(phi)*np.cos(phase[k])+z*np.sin(phi)*np.sin(phase[k]))/2, int(counts[k]), int(NumberOfCycles))
            # sometimes we still get nan ...
#            if probk[k] <= -40 or np.isnan(probk[k]):
            if probk[k] <= -40:
                probk[k] = -40
        else:
            probk[k] = -40
    lnprob = np.sum(probk)
    return lnprob

def lnbeta(p,f,N):
    # when p is veeery close to zero, we run into trouble, so we clip it at a 'reasonable' number
    p = np.clip(p,1./N**2,1-1./N**2)
    out= spec.gammaln(N+1)-spec.gammaln(1+f)-spec.gammaln(1+N-f)+f*np.log(p)+(N-f)*np.log(1-p)
    if np.isnan(out):
        print 'NaN detected'
        print (p,f,N)
    return out

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def get_errbars(prob, probx):
    maxind = np.argmax(prob)
    limit = .1
    psum = np.cumsum(prob)
    llim = find_nearest(psum, limit)
    ulim = find_nearest(psum, 1-limit)
    return probx[maxind], probx[llim], probx[ulim]

def bayesian_parity(dat, par_col=3):
    x = dat[:,0] / 360. 
    par = dat[:,par_col] 
    par = par - np.mean(par)
    y = (1+par) / 2.0
    a,b,c, p, px = fit_contrast_bayesian(x, y, NumberOfIons=1, mean_amp=.5, pm_mean_amp=0.5)
    maxlik, llim, ulim  = get_errbars(p,px)
    return maxlik, llim, ulim

if __name__ == "__main__":
    N = 100
    C = .05
    phi = np.pi/4.0
    x = np.linspace(0,1,20) 
    y = (1 + C * np.sin(x*  np.pi + phi)) / 2.0
    y1 = np.zeros_like(y)
    for i in xrange(y1.shape[0]):
        yy = np.random.multinomial(N,[y[i],1-y[i]],1) / float(N)
        y1[i] = yy[0][0]
    a,b,c, p, px = fit_contrast_bayesian(x, y1, NumberOfIons=1, mean_amp=.3, pm_mean_amp=0.3)
    maxlik, llim, ulim  = get_errbars(p,px)

