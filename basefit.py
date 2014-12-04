
from math import *
# 2011-01-11, TM: cluster only need to calculate, not plot ...
try:
    from pylab import *
except:
    pass
from scipy import *
import scipy.optimize as optimize
import numpy

"""basefit.py base class for fitting routines

This module contains FitBaseClass
"""

class FitBaseClass:
    plotobj = None
    default_take = "csingle"
    parameter_names = None
    """This is a base class for all fitting routines
    Look in fitdata.py how to use this base class
    """
    def __init__(self, data=None, take=None, use_ion=None, show_graph=True,
                 initial=None, quiet=False):
        self.data = None
        self.params = None
        self.fitfunc = None
        self.errfunc = lambda p, x, y: self.fitfunc(p,x) - y
        self.weighted_errfunc = lambda p, x, y, y_err: (self.fitfunc(p,x) - y) / y_err
        self.set_fitfunc()
        if not take:
            take = self.default_take
        self.take = take
        self.use_ion = use_ion
        self.show_graph = show_graph
        self.__get_data(data)
        self.fit_data(initial)
        if quiet == False:
            self.show_results()

    def set_fitfunc(self):
        raise RuntimeError("Fitting base class should not be called directly")

    def __get_data(self, my_data):
        """Reads and processes the datafile:
        If multiple ion data is found it tries to do a fit to each of the ions
        """
        if type(my_data) == numpy.ndarray:
            self.data = my_data
            return

    def get_initial_params(self, x, y, yerr):
        """Not implemented in the base class"""
        raise RuntimeError("This fitting routine does not support guessing of initial parameters")

    def fit_data(self, initial):
        """Fits the actual function
        If errors for the dataset are given it performs a weighted fit"""
        if self.data.shape[0] == 3:
            (x,y,yerr) = (self.data[:,0], self.data[:,1], self.data[:,2])
        else:
            (x,y) =  (self.data[:,0], self.data[:,1])
            yerr = None
        par_list = []
        try:
            self.NumberOfRows = y.shape[1]
        except IndexError:
            self.NumberOfRows = 1
        for index_rows in xrange(self.NumberOfRows):
            if self.NumberOfRows == 1:
                y_single = y
                yerr_single = yerr
            else:
                y_single = y[:,index_rows]
                yerr_single = yerr[:,index_rows]
            if not type(initial) == None:
                p0 = self.get_initial_params(x,y_single,yerr_single)
            else:
                p0 = initial
            if yerr == None:
                p1, covmat, infodict, mesg, ier = optimize.leastsq(self.errfunc, p0, args = (x, y_single),full_output=True)
            else:
                p1, covmat, infodict, mesg, ier = optimize.leastsq(self.weighted_errfunc, p0, args = (x, y_single, yerr_single),full_output=True)

            try:
                if len(p1) < 2:
                    raise RuntimeError()
            except:
                p1 = array([p1])
            chisq=sum(infodict["fvec"]*infodict["fvec"])
            dof=len(infodict["fvec"])-len(p1)
            try:
                
                sigma = sqrt(diag(covmat))*sqrt(chisq/dof)
            except:
                sigma = 0
            par_list.append((p1, sigma))


        self.param_dict = []
        for fit_index in xrange(self.NumberOfRows):
            my_dict = {}
            for param_index in xrange(len(self.parameter_names)):
                my_par = par_list[fit_index][0][param_index]
                try:
                    my_sigma = par_list[fit_index][1][param_index]
                except TypeError:
                    my_sigma = 100
                my_dict[self.parameter_names[param_index]] = (my_par, my_sigma)
            self.param_dict.append(my_dict)
        self.params = par_list

    def show_results(self):
        """Prints thex results"""
        for my_dict in self.param_dict:
            my_str = ""
            for key in my_dict:
                my_str += key + " - "
                my_str += str(my_dict[key][0]) + " - "
#            print my_str

        if self.show_graph:
            x = self.data[:,0]
            y = self.data[:,1]
            try:
                yerr = self.data[:,2]
            except:
                yerr = None
            
#            (x, y, yerr) = self.data
#            my_pdf = plotq.DefaultFrame
#            my_pdf.dataframe.new_axis("fit result")
            color_list = ['r','b','g','m','y','k','-r','-b','-g','-m','-y','-k']
            figure()
            for index_rows in xrange(self.NumberOfRows):
                my_color = color_list.pop(0)
                if self.NumberOfRows == 1:
                    y_single = y
                else:
                    y_single = y[:,index_rows]
#                hold(True)
                if index_rows == 0:
                    use_new_axis = True
                else:
                    use_new_axis = False
                if yerr != None:
                    yerr_single = yerr
                    errorbar(x, y_single, yerr_single, fmt='.'+my_color)
                else:
                    plot(x, y_single, '.'+my_color)

                x_new = linspace(min(x),max(x),500)
                plot(x_new, self.fitfunc(self.params[index_rows][0], x_new))
                show()

##
## basefit.py
## Login : <viellieb@odysseus>
## Started on  Sat Jan 30 21:36:29 2010
## $Id$
##
## Copyright (C) 2010
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

