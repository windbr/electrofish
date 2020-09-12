#!/usr/bin/env python3

__author__ = 'jpresern'

import sys
import numpy as np
from pyrelacs.DataClasses import load
from os import getcwd, walk
from utils import command_interpreter
from posixpath import join as ppjoin
from collections import OrderedDict
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from read_info_dat import read_info

def axes_def ():
    """

    :return fh: list of figure and axes handles
    """
    FHandles = plt.figure(figsize=(8.27, 6))
    ax1 = FHandles.add_axes([0.1, 0.38, 0.48, 0.53])

    return [FHandles, ax1]

def plot_spike_shapes(shape_dict, f_handle, meta, info, my_name):
    """
    Plots spike shapes superimposed, together with the appropriate mean
    :f_handle: list of figure and axes handles
    :param shape_dict: dictionary with spike shape for all selected repros and stimuli
    :return: figure handle, maybe
    """
    color_spread = np.linspace(0.8,0.4, len(shape_dict))
    cmapGrey = [ cm.Greys(x) for x in color_spread ]

    counter = 0
    for k, v in shape_dict.items():

        # for spike in range(shape_dict[k].shape[0]):
        # peak_voltage = np.max(shape_dict[k][spike,:])
        #     f_handle[1].plot(shape_dict[k][spike,:]-peak_voltage, color=cmapGrey[counter])

        avg_spike = np.mean(shape_dict[k],0)
        med_spike = np.median(shape_dict[k],0)

        peak_avg_voltage = np.max(avg_spike)
        peak_med_voltage = np.max(med_spike)

        q25spike = np.percentile(shape_dict[k],5,0)
        q75spike = np.percentile(shape_dict[k],95,0)
        std_spike = np.std(shape_dict[k],0)


        # f_handle[1].fill_between(range(0,med_spike.shape[0]), med_spike-peak_med_voltage-q25spike, med_spike-peak_med_voltage + q75spike, color=cmapGrey[counter], alpha=0.2)
        # f_handle[1].plot(med_spike-peak_med_voltage, color=cmapGrey[counter])

        f_handle[1].fill_between(range(0,avg_spike.shape[0]), avg_spike-peak_med_voltage-std_spike, med_spike-peak_med_voltage + std_spike, color=cmapGrey[counter], alpha=0.1)
        f_handle[1].plot(avg_spike-peak_avg_voltage, color=cmapGrey[counter])

        if "SyncPulse" in meta[k][0]["Status"].keys():
            f_handle[1].text(35, -5-counter*5, "".join([r'$g_{Vgate}$',' = ', meta[k][0]["Status"]["gvgate"].split('.')[0], ' ns']), fontsize=10, color=cmapGrey[counter])
            f_handle[1].text(60, -5-counter*5, "".join([r'$\tau_{Vgate}$',' = ', meta[k][0]["Status"]["vgatetau"].split('.')[0], ' ms']), fontsize=10, color=cmapGrey[counter])

        counter += 1

    f_handle[1].set_title(".".join(["_".join([info["Cell"]["Location"],expfolder,"spike_shape", my_name]), 'pdf']))

    f_handle[1].set_ylim(-40, 0)
    f_handle[1].set_xlabel("time [ms]")
    f_handle[1].set_ylabel("relative voltage [mV]")
    f_handle[1].set_xticklabels(['-1','0','1','2','3','4'])

    f_handle[0].savefig(ppjoin('../overviewSpikeShape/', ".".join(["_".join([info["Cell"]["Location"],expfolder,"spike_shape", my_name]), 'pdf'])), transparent=True)
    # f_handle[0].savefig(".".join(["Fig_spike_shape", 'pdf']), transparent=True)

def segment_extraction(fn1, fn2, ix, cur):
    """
    Detects spikes and extracts the segment around the spike, of desired length
    :param fn1: trace file name
    :param fn2: spikes file name
    :param index:
    :return spike_mat:2d matrix of spike shapes
    :return spike_time: 2d matrix of spike times
    """


    #   parse the Relacs file
    relacs_file1 = load(fn1)
    relacs_file2 = load(fn2)

    #   select the beloved section of Relacs file
    try:
        metas, _, spike_data = relacs_file2.select({"ReProIndex": int(ix), "I": cur})
        _, _, trace_data = relacs_file1.select({"ReProIndex": int(ix), "I": cur})
    except:
        return None

    shape_list = []

    for i in range(len(metas)):
        spikeses = spike_data[i][spike_data[i]>0]
        #   todo: make this going
        # spike_data[i][ 0 < spike_data[i] < np.float(metas[0]["Settings"]["Timing"]["duration"].split("ms")[0])]
        spike_wavelet = np.zeros([spikeses.shape[0], 100])
        for s in range(len(spikeses)):
            t_index = np.where(trace_data[i][:,0]==spikeses[s])
            spike_wavelet[s,:] = trace_data[i][t_index[0][0]-20:t_index[0][0]+80,1]

        shape_list.append(spike_wavelet)

    return np.vstack(shape_list), metas


if __name__ == "__main__":

    """
    Extracts the selected FI traces, detects spikes and plots them one over the other, to get average spike shape.
    Runs from command line in the current folder
    Input argument is a dictionary with the RePro indexes which you want to present in the same analysis as key.
    Stimulus current is used as a value
    example: report_spike_shape.py -d"{'16':'-0.14nA','34':'-0.14nA','40':'-0.14nA'}"
    """
    #   get the current working dir
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]

    #   get the command arguments
    ReProList = command_interpreter(sys.argv[1:])

    #   get the info.dat
    (_,_,filenames) = next(walk(ppjoin(wd, expfolder)))
    if "info.dat" in filenames:
        exp_info = dict(read_info(wd, expfolder)[0])
        print(exp_info)
    else:
        exp_info ={"Cell":{"Location":"UNLABELED"}}

    fn_trace = "ficurve-traces.dat"
    fn_spike = "ficurve-spikes.dat"

    #   compose list of RePro indices to be analysed
    ReProDict = ReProList['di']
    ReProIndex = sorted([k for k, v in ReProDict.items()])
    desired_shapes = OrderedDict()
    spike_shapes = OrderedDict()
    metka = OrderedDict()
    for r in ReProIndex:
        desired_shapes[r] = ReProDict[r]

    #   extract desired segment
    for ix, cur in desired_shapes.items():
        spike_shapes[ix], metka[ix] = segment_extraction(fn_trace, fn_spike, ix, cur)
    #   prepare name for the figure name
    name = '_'.join(map(str,ReProIndex))

    fig_handles = axes_def()

    plot_spike_shapes(spike_shapes, fig_handles, metka, exp_info, name)