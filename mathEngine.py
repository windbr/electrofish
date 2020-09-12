__author__ = 'jpresern'

import numpy as np
import peakutils
from math import log
from scipy.optimize import curve_fit

# from pylab import *

def func(x, a, c, d):

    """
    Defines several fitting mechanism
    """
    return a*np.exp(-c*x)+d

def exp_decay(timeList, voltageList, ySS=160):

    """
    This function finds a exponential decay
    : param tList : in sec
    : param vList : - voltage
    : param ySS : offset from which the exponential decay is computed
    : return parameters : a list of determined parameters (tau, amplitude ...)
    """

    parameters, _ = curve_fit(func, timeList, voltageList, p0=(1, 9, ySS), maxfev = 20000)
    # voltageFit = func(timeList, *parameters)
    return parameters

def exp_decay2(timeList, voltageList, p0=(220,90,160)):
    parameters, _ = curve_fit(func, timeList, voltageList, p0, maxfev = 20000)
    return parameters

def exp_fit(timeList, voltageList, ySS):

    """
    This function finds a
    :param tList : in sec
    :param vList : - voltage
    :param ySS : offset from which the exponential decay is computed
    :return amplitude : amplitude of the exp. fit
    :return tau : the time constant of the fit
    """

    bList = [log(max(y-ySS,1e-6)) for y in voltageList]
    b = np.matrix(bList).T
    rows = [ [1,t] for t in timeList]
    A = np.matrix(rows)
    #w = (pinv(A)*b)
    (w,residuals,rank,sing_vals) = np.linalg.lstsq(A,b)
    tau = -1.0/w[1,0]
    amplitude = np.exp(w[0,0])
    return (amplitude,tau)


def hist_peak_search(hist, bins):

    """
    Uses PeakUtil package to find the peak Locations and returns them as a value.
    It is used to determine the V_rest while the spikes are present
    : params hist : 1D array containing the hist midpoints and their amplitudes
    : params bins : 1D array containing the spread of data
    : returns peaks : peak locations (as value, not index)
    : retursn base : returns baseline, maybe it is more suitable for the V_rest evaluation

    """

    ix = peakutils.indexes(-hist, thres = 0.15/max(-hist), min_dist = 2)
    peaks = list(bins[list(ix)])

    return peaks

def baseline_search(voltage, pol_degree):

    """
    Uses PeakUtil package to find the baseline of the signal and returns it as a value.
    : params voltage : 1D array containing the input signal
    : params voltage : polynomial degree used in base line fitting
    : returns base : baseline of the signal
    """

    min = np.min(voltage)
    base = peakutils.baseline(voltage-min, pol_degree, 10000, 0.05)
    base = base+min
    return base