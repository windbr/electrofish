#!/usr/bin/env python3

__author__ = 'jpresern'

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys
from pyrelacs.DataClasses import load
from utils import read_stimuli_position
from os import getcwd
from utils import command_interpreter
from posixpath import join as ppjoin
from collections import OrderedDict

def resist_trace(path_to_folder, subfolder, raw_data_dict, metas):
    """
    Draws the resistance traces, one by one, from the raw trace
    :param raw_data_dict: Ordered dict with RePro indices as keys. Value is a two-element list containing iter of RePro
                            and np.array of index values.
    :param metas: metas pulled out from resistancetrace-exfit.dat
    :param subfolder: folder containing experimental data
    """

    #   few basic definitions
    dt = 0.05
    pauseDur = 50 + 10 # in miliseconds + 10 to get better figure

    #   how many RePro iterations?
    ReProIters = len(raw_data_dict)
    FHandles = plt.figure(figsize=(10, 3*ReProIters))

    #   trace figure dimensions
    trace_x = 0.1
    trace_width  = 0.3
    trace_height = 0.75/ReProIters
    trace_step = 1/ReProIters

    #   hist figure dimensions
    hist_x = 0.70
    hist_width = 0.25
    hist_height = 0.75/ReProIters
    hist_step = 1/ReProIters

    RePros = sorted([int(k) for k,v in raw_data_dict.items()])

    for i, k in enumerate(raw_data_dict):

        #   define and design grey color map
        color_spread = np.linspace(0.35,0.8, len(raw_data_dict[k]))
        cmapGrey = [ cm.Greys(x) for x in color_spread ]

        axL = FHandles.add_axes([trace_x,(i)*trace_step+(0.20/ReProIters), trace_width, trace_height])
        axR = FHandles.add_axes([hist_x,(i)*hist_step+(0.20/ReProIters), hist_width, hist_height])
        histData = np.zeros(1)
        # for j in range(raw_data_dict[k][1].shape[0]-1):
        for j in range(0, raw_data_dict[k].shape[0]-1):

            #   define time line
            timeLine = np.arange(-pauseDur, np.float(metas[int(k)]["Settings"]["Stimulus"]["duration"].split("m")[0]), dt)

            #   extract traces from raw file
            trace = np.fromfile("trace-1.raw", dtype=np.float32)[raw_data_dict[k][j,0] - np.int(pauseDur/dt) \
                        : np.int(np.float(metas[int(k)]["Settings"]["Stimulus"]["duration"].split("m")[0])/dt) + \
                    np.int(raw_data_dict[k][j,0])]
            #   plot traces from raw file
            axL.plot(timeLine, trace, color=cmapGrey[j])

            #   extract pre stimulus part from trace, compute mean and standard deviation, plot the lot.
            trace = trace[0:np.int(pauseDur/dt)]
            traceMean = np.mean(trace)


            #   offset everything to zero
            trace = trace - traceMean
            histData = np.concatenate((histData,trace))
            # histData = np.array(histData, dtype=np.float32)

            #   plot offseted hist data
            hist, bins = np.histogram(trace, 30, density=True)
            bincenters = 0.5*(bins[1:]+bins[:-1])
            axR.plot(bincenters,hist, color=cmapGrey[j])
            # axR.hist(trace,50, histtype="stepfilled", color=cmap[j])

        #   write text
        axL.text(0.40, 0.90, " ".join(["Experiment:", expfolder]), transform=axL.transAxes, fontsize=10)
        axL.text(0.40, 0.65, " ".join(['ReProIx:', str(k)]), transform=axL.transAxes, fontsize=10)
        axL.text(0.40, 0.55, " ".join(['DC current:', metas[int(k)]["Status"]["Current-1"]]), transform=axL.transAxes, fontsize=10)
        axL.text(0.40, 0.45, " ".join(['Amp. mode:', metas[int(k)]["Status"]["AmplifierMode"]]), transform=axL.transAxes, fontsize=10)

        # if "SyncPulse" in metas[int(k)]["Status"].keys():
        #     axL.text(0.05, 0.65, " ".join(['g:', metas[int(k)]["Status"]["g"]]), transform=axL.transAxes, fontsize=10)
        #     axL.text(0.05, 0.55, " ".join(['E:', metas[int(k)]["Status"]["E"]]), transform=axL.transAxes, fontsize=10)
        #     axL.text(0.05, 0.45, " ".join(['gVgate:', metas[int(k)]["Status"]["gvgate"]]), transform=axL.transAxes, fontsize=10)
        #     axL.text(0.05, 0.35, " ".join(['EVgate:', metas[int(k)]["Status"]["Evgate"]]), transform=axL.transAxes, fontsize=10)
        #     axL.text(0.05, 0.25, " ".join(['VgateTau:', metas[int(k)]["Status"]["vgatetau"]]), transform=axL.transAxes, fontsize=10)
        #     axL.text(0.05, 0.15, " ".join(['VgateMid:', metas[int(k)]["Status"]["vgatevmid"]]), transform=axL.transAxes, fontsize=10)
        #     axL.text(0.05, 0.05, " ".join(['VgateSlope:', metas[int(k)]["Status"]["vgateslope"]]), transform=axL.transAxes, fontsize=10)


        #   mean histogram
        histData = histData[1:]
        hist, bins = np.histogram(histData, 30, density=True)
        bincenters = 0.5*(bins[1:]+bins[:-1])
        axR.plot(bincenters, hist, color='k', linewidth=3)

        #   plot limits
        # axL.set_xlim(-pauseDur,np.float(metas[0]["Settings"]["Stimulus"]["duration"].split("m")[0]))
        axL.set_xlim(-pauseDur,20)
        axR.set_xlim(-2, 2)

        #   plot cosmetics
        axL.set_xlabel("Time [ms]")
        axL.set_ylabel("Voltage [mV]")
        axR.set_xlabel("Voltage [mV]")
        axR.set_ylabel("Count")
        axR.set_xticks([-2,-1,0,1,2])
        axL.axvline(-55, color='red')
        axL.axvline(-5, color='red')

        #   minimum statistics
        histSD = np.std(histData)

        #   write text
        axR.text(0.05, 0.8, " ".join(['signal std', str(np.round(histSD, 2))]), transform=axR.transAxes, fontsize=10)

    FHandles.savefig(ppjoin(".".join(["_".join([subfolder, "report_noise_traces", "_".join(map(str, RePros))]), 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join(["_".join([subfolder, "report_noise_traces", "_".join(map(str, RePros))]), 'svg'])), transparent=True)
    FHandles.savefig(ppjoin('../overviewNoise', ".".join(["_".join([subfolder,"report_noise_traces", "_".join(map(str, RePros))]), 'pdf'])), transparent=True)

if __name__ == "__main__":

    """
    Extracts the resistance traces from raw data, plots them in computes standard deviation of the signal on 50 ms of signal
    before the stimulus is delivered.
    Runs from command line in the experiment subfolder: report_noise_traces.py
    Input argument is a dictionary with the RePro indexes which you want to present in the same analysis as key. Stimulus iterations
    you want to exclude must be passed as a list into a dict value. If none, use [].
    example: report_noise_traces.py -d"{'2':[1,4,6], '4':[]}"
    """
    #   todo:   replace dict with ordered dict

    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]

    ReProList = command_interpreter(sys.argv[1:])

    filename = "membraneresistance-expfit.dat"
    relacs_file = load(filename)
    metas, _, _ = relacs_file.selectall()

    #   get all RePro executions
    ReProIxList = [str(metas[i]["ReProIndex"]) for i in range(len(metas))]

    #   convert metas into dict
    metas_dict = OrderedDict()
    for i in range(len(metas)):
        metas_dict[metas[i]["ReProIndex"]] = metas[i]

    #   read raw traces
    raw_data_dict = read_stimuli_position(wd, expfolder,'MembraneResistance', 2)

    #   replace the keys with the ReProIndices
    old_key_list = [k for k, v in raw_data_dict.items()]
    for i,key in enumerate(old_key_list):
        raw_data_dict[ReProIxList[i]] = raw_data_dict.pop(key)

    #   select RePro demanded by command
    demanded_RePros = sorted([int(k) for k, v in ReProList["od"].items()])
    selected_data_dict = OrderedDict()
    for k in demanded_RePros:
        selected_data_dict[str(k)] = raw_data_dict[str(k)]

    #   select raw traces according to the command, default is to exclude the last
    for k,v in selected_data_dict.items():
        exclude = ReProList["od"][k]
        if exclude:
            selected_data_dict[k]= np.delete(selected_data_dict[k], exclude, axis=0)
    #   draw and measure
    resist_trace(wd, expfolder, selected_data_dict, metas_dict)