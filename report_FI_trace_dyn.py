#!/usr/bin/env python3

__author__ = 'janez'

import sys, getopt, ast, types
import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from os import walk, getcwd
from collections import OrderedDict, defaultdict
from IPython import embed
from utils import command_interpreter
from FI import fi_spikes_sort, fi_instant_freq, fi_inst_freq_continuous

# return (arg_dict["li"])

def fi_voltage_plot(ReProList,iterations):

    #   defines relative sizes of individual raster plots
    rug_bot = 0.1
    rug_width  = 0.4
    rug_height = 0.67/(len(ReProList))
    rug_step = 1/(len(ReProList))

    #   defines the color map
    cmapBlue = ["Blue", "DarkTurquoise","CadetBlue", "DeepSkyBlue", "CornFlowerBlue", "DodgerBlue", "LightSkyBlue", "LightSteelBlue", "DarkCyan"]

    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]

    #   define the input data file name
    filename = "ficurve-traces.dat"

    #   get the file and load it
    relacs_file = load(filename)

    #   define the absolute size of the raster & return plot
    FHandles = plt.figure(figsize=(14, 3*len(ReProList) + 2))

    #   counter
    counter = 0

    #   define OrderedDict for freq_continuous between repros.
    freq = OrderedDict()

    #   define g_list
    g_list = []
    tau_list = []

    for k, v in ReProList.items():
        metas, _, datas = relacs_file.select({"ReProIndex": int(k), "I": v})

        trace = datas[iterations[0]][:,1]
        timeLine = datas[iterations[0]][:,0]

        #   add axes
        axarr = FHandles.add_axes([rug_bot,counter*rug_step+(0.25/(len(ReProList)+1)), rug_width, rug_height])
        axarr.plot(timeLine, trace, '|', color=cmap[counter])

        axarr.set_xlabel('Time [ms]')
        axarr.set_ylabel('Voltage [ms]')

        if "SyncPulse" in metas[0]["Status"].keys():
            g_list.append(metas[0]["Status"]["gvgate"])
            tau_list.append(metas[0]["Status"]["vgatetau"])
            #   draws R against g
            axarr.text(0.02, 0.15, " ".join(['gvgate = ', metas[0]["Status"]["gvgate"]] ), color=cmap[counter], transform=axarr.transAxes, fontsize=10)
            axarr.text(0.8, 0.15, " ".join(['Vtau = ', metas[0]["Status"]["vgatetau"]] ), color=cmap[counter], transform=axarr.transAxes, fontsize=10)
        #   grid
        axarr.grid(b=True)

        axarr.set_xlim(0,400)

        #   increase counter
        counter = counter + 1

    ##   prepare file names
    #   extract ReProIx, sort them and merge them
    keys  =[i for i in ReProList]
    keys = sorted(list(map(int, keys)))
    keys = '.'.join(map(str, keys))
    #   extract current values
    current = list(ReProList.values())[0]
    #   join keys, currents and iterations
    name ='.'.join([keys, current, 'iter', str(iterations[0])])

    #   Save figures
    FHandles.savefig(ppjoin(".".join([expfolder, "dyn_FI_trace", name, 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join([expfolder, "dyn_FI_trace", name, 'svg'])), transparent=True)

if __name__ == "__main__":

    """
    Shows the raw traces of a single FI repro iteration (selected) at certain stimulus current.
    Runs from command line in the experiment subfolder: report_FI_trace_dyn.py
    Another input argument is a list with a *single* iteration number
    Input argument is a dictionary with a RePro indexes as a key and current as a value
    example: report_FI_trace_dyn.py -l [3] -d "{'291': '-0.25nA', '293': '-0.25nA'}"
    """
    #   get the current working dir, split it into the dir and the path to it
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]
    ReProList = command_interpreter(sys.argv[1:])
    print(ReProList["di"])
    fi_voltage_plot (ReProList["di"], ReProList["li"])


