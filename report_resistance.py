#!/usr/bin/env python3

__author__ = 'janez'

import sys, getopt, ast, types
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from mathEngine import hist_peak_search, exp_decay
from os import walk, getcwd
from collections import OrderedDict, defaultdict
from IPython import embed
from utils import command_interpreter

def aling_to_zero(trace_orig):
    """
    Takes the voltage trace before stimulus, computes the Vrest and resets the whole trace to zero. Returns the new traces.
    :param trace_orig: np.array containing voltage trace
    :param
    :return trace_norm: np.array normalized voltage trace
    """
    #   histogram calculation for the Vrest distribution
    hist, bins = np.histogram(trace_orig[0, trace_orig[0,:,0]<0, 1], 50, density=True)
    bincenters = 0.5*(bins[1:]+bins[:-1])

    #   find a V_rest
    memPotCalc = hist_peak_search(-hist, bincenters)

    #   if V_rest fails, try to compute median
    if memPotCalc == []:
        memPotCalc = [0,0]
        memPotCalc[0] = np.median(trace_orig[0, trace_orig[0,:,0]<0, 1])

    #   subtracts calculated V_rest from the whole trace
    trace_norm = trace_orig
    trace_norm[0,:,1] = trace_orig[0,:,1] - memPotCalc[0]

    #   returns corrected trace
    return trace_norm

def computes_resistance(trace, amplitude, duration):
    """
    Takes current and voltage traces, extracts the steady-state part of both, computes median and using R=U/I computes the resistance
    :param trace: np array containing signal traces
    :param duration: number containing duration of the signal
    """

    # current = np.median(trace[0, (trace[0,:,0]>duration*0.7) & (trace[0,:,0]<duration), 3])
    current = abs(amplitude)
    voltage = np.median(trace[0, (trace[0,:,0]>duration*0.7) & (trace[0,:,0]<duration), 1])
    R = voltage/current
    return abs(np.round(R))


def resist_plot(ReProIx, wd, expfolder, norm=False):
    """
    Takes resistances from selected RePro iterations and plots them one over the other
    :param ReProIx: list of indexes you want to plot
    :param norm: whether or not you want take difference in offsets into account
    """

    #   define file name
    filename = "membraneresistance-trace.dat"

    #   load data
    relacs_file=load(filename)

    # #   four panel figure
    # FHandles = plt.figure(figsize=(10, 10))
    # axarr_orig = FHandles.add_axes([0.1, 0.7, 0.85, 0.25])
    # axarr = FHandles.add_axes([0.10, 0.375, 0.85, 0.25])
    # axarrR = FHandles.add_axes([0.60, 0.05, 0.35, 0.25])
    # axarrTau = FHandles.add_axes([0.10, 0.05, 0.35, 0.25])

    #   three panel figure
    FHandles = plt.figure(figsize=(10, 8))
    axarr = FHandles.add_axes([0.10, 0.55, 0.85, 0.4])
    axarrR = FHandles.add_axes([0.60, 0.07, 0.35, 0.4])
    axarrTau = FHandles.add_axes([0.10, 0.07, 0.35, 0.4])

    #   define the colormap
    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
             "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]
    # cmapBlue = ["DarkTurquoise","DarkCyan", "CadetBlue", "DeepSkyBlue", "CornFlowerBlue", "DodgerBlue", "LightSkyBlue", "LightSteelBlue"]

    #   counter
    counter = 0

    #   define trace list
    V_trace_list = []

    #   define resistance list
    R_list = []

    #   define g list
    g_list = []

    #   define tau list
    tau_list = []

    for ix in ReProIx:

        #   if relacs file is empty or too short due to aborted RePro presentation
        #   "try:" provides normal termination of the loop instead of error
        try:
            metas, _, datas = relacs_file.select({"ReProIndex": ix})

        except:
            return None
        #   convert into np array
        trace = np.array(datas)

        #   extract stimulus duration
        durStim = float(metas[0]["Settings"]["Stimulus"]["duration"].split('m')[0])

        #   plot original, non-normalized trace
        # axarr_orig.plot(trace[0,:,0], trace[0,:,1], color=cmap[counter])

        #   normalize to the 0 mV
        if norm==True:
            trace = aling_to_zero(trace)

        #   calculate the resistance from voltage drop in the steady-state part
        R = computes_resistance(trace, float(metas[0]["Settings"]["Stimulus"]["amplitude"].split('n')[0]), float(metas[0]["Settings"]["Stimulus"]["duration"].split('m')[0]))
        R_list.append(R)

        #   fits exponential
        parameters = exp_decay(trace[0, (durStim > trace[0,:,0]) & (0<trace[0,:,0]), 0],
                               trace[0, (durStim > trace[0,:,0]) & (0<trace[0,:,0]), 1], 0)

        fitVoltage = parameters[0]*np.exp(-parameters[1]*trace[0, (durStim > trace[0,:,0]) & (0<trace[0,:,0]), 0])\
                     +parameters[2]

        #   append tau_list
        tau_list.append(1/parameters[1])

        #   print RePro iteration
        print("ReProIx", ix, "Iterations", len(metas))

        #   append trace list
        V_trace_list.append(trace[0,:,1])

        #   draw within loop
        axarr.plot(trace[0,:,0], trace[0,:,1], color=cmap[counter])

        #   extracts the g & plots
        if "SyncPulse" in metas[0]["Status"].keys():
            g_list.append(float(metas[0]["Status"]["g"].split('n')[0]))
            #   draws R against g
            axarrR.plot(float(metas[0]["Status"]["g"].split('n')[0]), R, 'o', color=cmap[counter])
            axarrTau.plot(float(metas[0]["Status"]["g"].split('n')[0]), 1/parameters[1], 'o', color=cmap[counter])
            axarr.plot(trace[0, (durStim > trace[0,:,0]) & (0<trace[0,:,0]), 0], fitVoltage, '--', color='k')
            axarr.text(0.05, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["g"]] ), color=cmap[counter], transform=axarr.transAxes, fontsize=10)

        #   counter increase
        counter = counter+1

    #   draws plot, flips and transforms data holders
    # time=trace[0,:,0]
    # V_trace_list = np.array(V_trace_list)
    # V_trace_list = V_trace_list.T
    # axarr.plot(time, V_trace_list)
    # axarr.fill_between(trace[0,:,0], trace[0,:,1]+trace[0,:,2], trace[0,:,1]-trace[0,:,2], color=cmapBlue[counter], alpha=0.2)

    #   resistance
    print(R_list)

    #   writes x labels
    axarr.set_xlabel('Time [ms]')

    #   writes y labels
    axarr.set_ylabel('Relative voltage [mV]')

    #   grid
    axarr.grid(b=True)

    #   plot comments and other annotations
    axarr.text(0.25, 0.90, " ".join(['Apteronotus leptorhynchus', expfolder]), transform=axarr.transAxes, fontsize = 10)
    axarr.text(0.65, 0.15, " ".join(['Stimulus amplitude:', metas[0]["Settings"]["Stimulus"]["amplitude"]] ), transform=axarr.transAxes, fontsize=10)
    axarr.text(0.65, 0.05, " ".join(['Stimulus duration:', metas[0]["Settings"]["Stimulus"]["duration"]] ), transform=axarr.transAxes, fontsize=10)
        # axarrR.text(0.25, 0.15, " ".join(['Resistance:', comments4["Rss"]]), transform=axarrR.transAxes, fontsize = 10, color = 'r')

    #   draws R against g
    if "SyncPulse" in metas[0]["Status"].keys():
        # axarrR.plot(g_list, R_list, 'o')
        axarrR.set_ylabel('Membrane resistance [MOhm]')
        axarrR.set_xlabel('g leak [nS]')
        axarrTau.set_ylabel('Tau [ms]')
        axarrTau.set_xlabel('g leak [ns]')

    #   sets lims
    axarrR.set_xlim(min(g_list)*1.1, max(g_list)*1.1)
    axarrTau.set_xlim(min(g_list)*1.1, max(g_list)*1.1)
    axarrR.set_ylim(0, max(R_list)*1.1)
    axarrTau.set_ylim(-1, max(tau_list)*1.1)

    #   Save figures
    FHandles.savefig(ppjoin(".".join([expfolder, "resistance", 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join([expfolder, "resistance", 'svg'])), transparent=True)

if __name__ == "__main__":

    """
    Runs from command line in the experiment subfolder: report_resist.py
    Input argument is a list with the RePro indexes which you want to present in the same analysis.
    example: report_resist.py -l [15,17,21,23]
    """
    #   get the current working dir, split it into the dir and the path to it
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]
    ReProList = command_interpreter(sys.argv[1:])
    resist_plot(ReProList["li"], wd=wd, expfolder=expfolder, norm=True)