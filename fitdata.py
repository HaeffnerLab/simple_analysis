"""
fitdata.py

This file contains advanced fitting functions which are based on the FitBaseClass

Available functions:

SineFit
RabiFit - initial is only the Rabi time
PeakFit
DipFit
ExpFit
BlueFit
GaussianFit
MultiPeakFit        
RabiDecayFit

Syntax:
-------

FunctionName(data=dataobj/timestring/array, take=None, use_ion=None,
             show_graph=False, initial=None)

Params:
-------

take - one of "qc", "csingle", "cprb"
use_ion - If this is set then just this ion is used for the csingle data.
show_graph - Show a graph with the fitted values - NOT Implemented yet
initial - Initial parameters. If omitted the fit functions tries to guess the initial parameters


Example:
--------
Load Datapath and fit a sine function to the dataset qc1815
>>> aus(20100118)
>>> f1 = fit.SineFit('1815')

If you want to see the result you should type:

>>> f1 = fit.SineFit('1815', show_graph=True)

To use the results:

A tuple with amplitude and its error from the first ion:
>>> f1.param_dict[0]['Amplitude']

The TimeConstant of the second ion:
>>> f1.param_dict[1]['TimeConstant'][0]

Get all parameters for ion n in an array:
>>> f1.params[n][0]

Get all errors for ion n in an array:
>>> f1.params[n][1]

----------------------------
Depreciated - DO NOT USE
Some fits may also be called with the function

fit.FitRoutines
"""

from math import *
# 2011-01-11, TM: cluster only need to calculate, not plot ...
from pylab import *
from scipy import *
import scipy.optimize as optimize

from basefit import FitBaseClass

class SineFit(FitBaseClass):
    """Fits a sine to the data
    """
    default_take = "csingle"
    parameter_names = ['Amplitude', 'TimeConstant', 'Phase', 'Offset']

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: p[0]/2*cos(2*pi/p[1]*x+p[2]) + p[3]
        
    def get_initial_params(self, x, y, yerr):
        """Guess initial parameters with a FFT"""
        spec = fft(y)
        freq_spec = fftfreq(len(y),d=x[1]-x[0])
        estimated_freq = max(freq_spec[find(max(abs(spec[:,nonzero(freq_spec>0)[0]]))==abs(spec))])
        estimated_time = 1./estimated_freq
        p0 = array([max(y)-min(y), estimated_time, 0, (max(y)+min(y))/2])
        return p0


class SineFitPosAmpl(SineFit):
    """Only positive amplitudes allowed """
    def set_fitfunc(self):
        self.fitfunc = lambda p, x: abs(p[0])/2*cos(2*pi/p[1]*x+p[2]) + p[3]

class RabiFit(FitBaseClass):
    """RabiFit is just a sine fit with different guesses for the amplitude and offset
    the optional parameter initial gives an initial parameter for the rabitime
    """
    default_take = "csingle"
    parameter_names = ['TimeConstant', 'Phase', 'Offset']

#   def fit_data(self,initial):
#       """We rewrite the initial variable"""
#       if initial != None:
#           initial = array([initial,0,.5])
#
#       SineFit.fit_data(self,initial)

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: 1.0/2*cos(2*pi/p[0]*x+p[1]) + p[2]


    def get_initial_params(self, x, y, yerr):
        """Set initial guess for amplitude and offset"""
        spec = fft(y)
        freq_spec = fftfreq(len(y),d=x[1]-x[0])
        estimated_freq = max(freq_spec[find(max(abs(spec[:,nonzero(freq_spec>0)[0]]))==abs(spec))])
        estimated_time = 1./estimated_freq
        p0 = array([estimated_time, 0, .5])
        return p0

class RabiDecayFit(FitBaseClass):
    """RabiDecayFit fits an exponentially decaying Rabi flop
    """
    default_take = "csingle"
    parameter_names = ['RabiTimeConstant', 'Phase', 'Offset', 'DecayTimeConstant']

#   def fit_data(self,initial):
#       """We rewrite the initial variable"""
#       if initial != None:
#           initial = array([initial,0,.5])
#
#       SineFit.fit_data(self,initial)

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: exp(-x/abs(p[3])) * 1.0/2*cos(2*pi/p[0]*x+p[1]) + p[2]


    def get_initial_params(self, x, y, yerr):
        """Set initial guess for amplitude and offset"""
        spec = fft(y)
        freq_spec = fftfreq(len(y),d=x[1]-x[0])
        estimated_freq = max(freq_spec[find(max(abs(spec[:,nonzero(freq_spec>0)[0]]))==abs(spec))])
        estimated_time = 1./estimated_freq
        p0 = array([estimated_time, 0, .5, 1000.])
        return p0


class LinFit(FitBaseClass):
    """Dip fit fits an inverse Lorentzian"""
    default_take = "qc"
    parameter_names = ['Slope',"Offset"]

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: p[0] * x + p[1]

    def get_initial_params(self, x, y, yerr):
        """Find the position of the minimum value"""
        p = [0,0]
        p[1] = min(y)
        p[0] = (max(y)-min(y))/(max(x)-min(x))
        return array(p)


class DipFit(FitBaseClass):
    """Dip fit fits an inverse Lorentzian"""
    default_take = "qc"
    parameter_names = ['Center', 'Width', 'Asympt. Maximum', 'Minimum']

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: (p[2] - (p[2]-p[3])*(p[1]/2)**2 / ((p[0] - x)**2 + p[1]/2) **2)

    def get_initial_params(self, x, y, yerr):
        """Find the position of the minimum value"""
        estimated_max = max(y)
        estimated_min = min(y)
        y1 = map(int, y *1000)
        estimated_position = x[ y1.index(min(y1)) ]
        estimated_width = (max(x) - min(x)) / 20.0
        p0 = array([estimated_position, estimated_width, estimated_max, estimated_min])
        return p0


class PeakFit(FitBaseClass):
    """Fits a Lorentzian to the data"""
    default_take = "qc"
    parameter_names = ["Center", "Width", "Amplitude"]

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: (p[2]) * (p[1]/2)**2 / ( (p[0] -x )**2 + (p[1]/2)**2 )

    def get_initial_params(self, x, y, yerr):
        """Find the position of the maximum value"""
        estimated_height = max(y)
        y1 = map(int, y *1000)
        estimated_position = x[ y1.index(max(y1)) ]
        estimated_width = (max(x) - min(x)) / 20.0
        p0 = array([estimated_position, estimated_width, estimated_height])
        return p0


class GaussianFit(FitBaseClass):
    """Fits a Gaussian to the data"""
    default_take = "qc"
    parameter_names = ["Center", "Width", "Amplitude"]

    def set_fitfunc(self):
        self.fitfunc = lambda p, x: p[2] * exp(- (x - p[0]) ** 2. / (2 * p[1] ** 2))

    def get_initial_params(self, x, y, yerr):
        """Find the position of the maximum value"""
        estimated_height = max(y)
        y1 = map(int, y *1000)
        estimated_position = x[ y1.index(max(y1)) ]
        estimated_width = (max(x) - min(x)) / 20.0
        p0 = array([estimated_position, estimated_width, estimated_height])
        return p0


class MultiPeakFit(GaussianFit):
    '''Fits multiple peaks to the data.'''
    default_take = "csingle"
    

def fitfunc_blue(p,x):
    y = zeros(x.shape[0])
    for i in xrange(10):
#        y += abs(p[i+1])/2.0 * (1-cos(2*pi/p[0]*sqrt(i+1) * x))
#        y += abs(p[1])**i/(p[1]+1)**(i+1) /2.0 * (1-cos(2*pi/p[0]*sqrt(i+1) * x))
        y += abs(p[1])**i/(p[1]+1)**(i+1) /2.0 * (1-cos(2*pi/90*sqrt(i+1) * x))

    return y


class BlueFit(FitBaseClass):
    """Fits a thermal occupation to the blue sideband"""
    default_take = "qc"
    parameter_names = ["Time",'MeanPhonon']
    max_phonons = 10


    def set_fitfunc(self):
        
#        if self.parameter_names[-1][0] != 'p':
#            for i in xrange(self.max_phonons-1):
#                self.parameter_names.append('phonon'+str(i))
        self.fitfunc = fitfunc_blue

    def get_initial_params(self, x, y, yerr):
        """Find the position of the maximum value"""
#        p0 = zeros(self.max_phonons + 1)
        p0 = zeros(2)
        p0[0] = 100
        p0[1] = .1
        return p0



class ExpFit(FitBaseClass):
    """Fits an exponential decay
    Parameters
      use_offset: if set to True then a decay with a constant offset is fitted
      amplitude: if not None then this is the fixed amplitude Does not support offset
      fixed offset: is used if used_offset is false
    """
    default_take = "qc"

    def __init__(self, data=None, take=None, use_ion=None,
                 show_graph=True, initial=None, use_offset = False , 
                 fixed_offset= 0.0, amplitude=None):
        self.offset = use_offset
        self.fixed_offset = fixed_offset
        self.amplitude = amplitude
        if amplitude != None:
            self.parameter_names = ["TimeConstant"]
        else:
            if use_offset:
                self.parameter_names = ["TimeConstant", "Amplitude", "Offset"]
            else:
                self.parameter_names = ["TimeConstant", "Amplitude"]
        FitBaseClass.__init__(self, data=data, take=take, use_ion=use_ion,
                 show_graph=show_graph, initial=initial)

    def set_fitfunc(self):
        """Set the function with or without offset"""
        if self.amplitude != None:
#            print self.amplitude
            self.fitfunc = lambda p, x: (self.amplitude * exp(-x * p[0]))
        else:            
            if self.offset:
                self.fitfunc = lambda p, x: (p[1] * exp(-x * p[0]) + p[2])
            else:
                self.fitfunc = lambda p, x: (p[1] * exp(-x * p[0]) + self.fixed_offset)

    def get_initial_params(self, x, y, yerr):
        """Get a guess on the parameters"""
        ampl = y[0]
        offset = 0
        tau = log(y[-1] / float(y[0])) / (x[-1] - x[0])
        if self.amplitude != None:
            p0 = array([tau])
        else:
            if self.offset:
                p0 = array([tau, ampl, offset])
            else:
                p0 = array([tau, ampl])
        return p0

class FitRoutines:
    def __init__(self, data, type=None, take="csingle", use_ion=None, show_graph=False,
                 initial=None):
        self.type_dict={}
        self.type_dict['sine'] = SineFit
        self.type_dict['peakfit'] = PeakFit
        self.fitobj = self.type_dict[type](data, take=take, use_ion=use_ion,
                                           show_graph=show_graph, initial=initial)
        self.params = self.fitobj.params
        params = self.params[0][0]
        sigma = self.params[0][1]
        self.param_dict = self.fitobj.param_dict
#        self.param_dict['Amplitude'] = (params[0],sigma[0])
#        self.param_dict['TimeConstant'].append((params[1],sigma[1]))
#        self.param_dict['Phase'].append((params[2],sigma[2]))
#        self.param_dict['Offset'].append((params[3],sigma[3]))


##
## fitdata.py
## Login : <viellieb@ohm>
## Started on  Sat Aug  2 11:33:56 2008 Philipp Schindler
## $Id$
##
## Copyright (C) 2008 Philipp Schindler
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##


# fitdata.py ends here
