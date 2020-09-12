#!/usr/bin/env python3

__author__ = 'jpresern'

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from os import walk, getcwd
from collections import OrderedDict
from utils import command_interpreter
from FI import fi_spikes_sort, fi_instant_freq, fi_inst_freq_continuous
from read_info_dat import read_info

# return (arg_dict["li"])

def fi_raster_return(ReProList_Unordered, info):

    #   Order!
    ReProList = OrderedDict(sorted(ReProList_Unordered.items(), key=lambda x:x[0], reverse=False))

    #   defines relative sizes of individual raster plots
    rug_bot = 0.1
    # rug_width  = 0.48
    rug_width  = 0.22
    rug_height = 0.67/(len(ReProList) + 1)
    rug_step = 1/(len(ReProList) + 1)

    #   defines relative size of individual return plots
    return_bot = 0.67
    return_width = 0.11
    return_height = 0.67/(len(ReProList) + 1)
    return_step = 1/(len(ReProList) + 1)

    #   defines relative size of individual histograms
    hist_bot =  0.85
    hist_width = 0.11
    hist_height = 0.67/(len(ReProList) + 1)
    hist_step = 1/(len(ReProList) + 1)

    #   defines the color map
    cmapBlue = ["Blue", "DarkTurquoise","CadetBlue", "DeepSkyBlue", "CornFlowerBlue", "DodgerBlue", "LightSkyBlue", "LightSteelBlue", "DarkCyan"]

    cmap = ["DarkOrange", "Blue", "LimeGreen", "Gold", "Plum", "DarkRed", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]


    #   define and design grey color map
    color_spread = np.linspace(0.35,0.9, len(ReProList))
    cmapGrey = [ cm.Greys(x) for x in color_spread ]

    #   define the input data file name
    filename = "ficurve-spikes.dat"

    #   get the file and load it
    relacs_file = load(filename)

    #   define the absolute size of the raster & return plot
    FHandles = plt.figure(figsize=(18, 3*len(ReProList) + 3))
    # FHandles = plt.figure(figsize=(18, 3*len(ReProList) + 3))

    #   counter
    counter = 0

    #   define OrderedDict for freq_continuous between repros.
    freq = OrderedDict()

    #   define g_list
    g_list = []
    tau_list = []

    for k, v in ReProList.items():
        metas, _, datas = relacs_file.select({"ReProIndex": int(k), "I": v})

        #   sorts spikes
        spikeDict, stim_amps = fi_spikes_sort(metas, datas)

        #   computes instantaneous spike frequencies based on spikes
        freqDict, timeDict = fi_instant_freq(spikeDict) #, metas[counter[i]:counter[i+1]])
        freq_continuous, _, timeLine, freqSD_continuous = fi_inst_freq_continuous(spikeDict, freqDict, metas)

        #   add axes
        axarr = FHandles.add_axes([rug_bot,counter*rug_step+(0.25/(len(ReProList)+1)), rug_width, rug_height])
        axarr2 = axarr.twinx()
        #   plot
        axarr2.plot(timeLine, freq_continuous[v], color = cmapGrey[counter])
        for j in range(len(spikeDict[v])):
            axarr.plot(spikeDict[v][j], np.ones(len(spikeDict[v][j]))+j, '|', color=cmap[j], linewidth=2, mew=1)

        axarr.set_xlabel('Time [ms]')
        axarr.set_ylabel('Iter')
        axarr2.set_ylabel('Frequency [Hz]', color=cmap[counter])

        #   insert computed freq into OrderedDict
        freq[k] = freq_continuous[v]

        if "SyncPulse" in metas[0]["Status"].keys():
            g_list.append(metas[0]["Status"]["gvgate"].split(".")[0])
            tau_list.append(metas[0]["Status"]["vgatetau"].split(".")[0])
            #   draws R against g r'$\alpha_i > \beta_i$'
            axarr.text(0.02, 0.8, " ".join([r'$g_{Vgate}$',' = ', metas[0]["Status"]["gvgate"].split('.')[0],' nS' ]), color=cmapGrey[counter], transform=axarr.transAxes, fontsize=12)
            axarr.text(0.7, 0.8, " ".join([r'$\tau_{Vgate}$',' = ', metas[0]["Status"]["vgatetau"].split('.')[0],' ms' ] ), color=cmapGrey[counter], transform=axarr.transAxes, fontsize=12)
        #   grid
        axarr.grid(b=True)

        #   set y lim
        axarr.set_ylim(0, len(spikeDict[v])+1)
        axarr.set_xlim(-100,500)

        #   drawing the second plot, containing time intervals
        axarrM = FHandles.add_axes([return_bot,counter*hist_step+(0.25/(len(ReProList)+1)), return_width, return_height])

        allIntervalsZero = np.array([0])
        for gnj in range(len(timeDict[v])):
            allIntervals = np.append(allIntervalsZero, timeDict[v][gnj])

            for m in range(len(allIntervals)-1):
               # axarrM.plot(allIntervals[m],allIntervals[m+1], "o", color= cmap[gnj], ms=5)
               #  axarrM.plot(allIntervals[m],allIntervals[m+1]-allIntervals[m], "o", color= cmap[gnj], ms=5)
                axarrM.plot(allIntervals[m],allIntervals[m+1]-allIntervals[m], "o", color = cmapGrey[counter], ms=5)

        # axarrR.set_xlim(0, np.percentile(allIntervals,97))
        # axarrR.set_ylim(0, np.percentile(allIntervals,97))

        axarrM.set_xlim(0, 12)
        axarrM.set_ylim(-6, 6)
        # axarrM.set_xlim(0, np.percentile(allIntervals,97))
        # axarrM.set_ylim(0, np.percentile(allIntervals,97))

        axarrM.set_xlabel('n ISI [ms]')
        axarrM.set_ylabel('n+1 ISI - n ISI [ms]')

        #   drawing the third plot, containing time interval distribution
        axarrR = FHandles.add_axes([hist_bot,counter*hist_step+(0.25/(len(ReProList)+1)), hist_width, hist_height])

        #   put all insterspike intervals into a single numpy array
        allIntervals = []

        for znj in range(len(timeDict[v])):
            IntervalsN = np.append(allIntervalsZero, timeDict[v][znj])
            # for ml in range((IntervalsN).shape[0]-1):
            #     allIntervals.append((IntervalsN[ml+1]-IntervalsN[ml]))

            IntervalsN1 = np.append(timeDict[v][znj], allIntervalsZero)
            allIntervals.extend((IntervalsN1-IntervalsN)[1:-2])

        allIntervals = np.array(allIntervals)
        # allIntervals = allIntervals[1:-1]
        #   histogram calculation for the ISI distribution
        if allIntervals.shape[0] > 1:
            # allIntervals = allIntervals[allIntervals<np.percentile(allIntervals,95)]
            if allIntervals.shape[0] > 1:
                axarrR.hist (allIntervals, bins = 20, normed = False, color= cmapGrey[counter], histtype='stepfilled', orientation='horizontal')

        #   add labels
        axarrR.set_ylabel('n+1 ISI - n ISI [ms]')
        axarrR.set_xlabel('Count')

        #   set x lim
        # axarrR.set_xticks([0, 5, 10, 15])
        # axarrR.set_xticklabels([0, 5, 10, 15])
        # axarrR.set_xlim(0, np.percentile(allIntervals,97))
        axarrR.set_ylim(-6, 6)

        #   increase counter
        counter = counter + 1

    #   add top plot to compare instantaneous frequencies
    axarr = FHandles.add_axes([rug_bot,counter*rug_step+(0.25/(len(ReProList)+1)), rug_width, rug_height])
    for l, w in enumerate(freq):
        axarr.plot(timeLine, freq[w], color = cmapGrey[l])

    #   writes x labels
    axarr.set_xlabel('Time [ms]')

    #   writes y labels
    axarr.set_ylabel('Frequency [Hz]')

    #   grid
    # axarr.grid(b=True)

    axarr.set_xlim([-50,550])

    #   plot comments and other annotations
    axarr.text(0.25, 0.90, " ".join(['Apteronotus leptorhynchus', expfolder]), transform=axarr.transAxes, fontsize = 10)
    axarr.text(0.55, 0.75, " ".join(['Stimulus amplitude:', metas[0]["I"]] ), transform=axarr.transAxes, fontsize=10)

    #   extracts the g & plots
    if "SyncPulse" in metas[0]["Status"].keys():
        #   draws R against g
        for i in range(len(g_list)):
            axarr.text(0.75, 0.60-0.1*i, " ".join([r'$\tau_{Vgate}$',' = ', tau_list[i], ' ms']), color=cmapGrey[i], transform=axarr.transAxes, fontsize=10)
            axarr.text(0.02, 0.60-0.1*i, " ".join([r'$g_{Vgate}$',' = ', g_list[i], ' nS']), color=cmapGrey[i], transform=axarr.transAxes, fontsize=10)

    #   extract ReProIx, sort them and merge them
    keys  =[i for i in ReProList]
    keys = sorted(list(map(int, keys)))
    keys = '.'.join(map(str, keys))
    #   extract current values
    current = list(ReProList.values())[0]
    #   join keys, currents and iterations
    name ='_'.join([keys, current])

    #   write title on the figure
    FHandles.suptitle("".join([info["Cell"]["Location"],':', expfolder, "FI_freq_dyn", name]), fontsize=12)

    #   Save figures
    FHandles.savefig(ppjoin(".".join([expfolder, "FI_freq", name, 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join([expfolder, "FI_freq", name, 'svg'])), transparent=True)
    #   dump .pdf to selected folder
    FHandles.savefig(ppjoin('../overviewFI/', ".".join(["_".join([info["Cell"]["Location"],expfolder,"FI_freq_dyn", name]), 'pdf'])), transparent=True)

if __name__ == "__main__":

    """
    Plots spike frequency comparison between between selected RePro executions and stimuli.
    Runs from command line in the experiment subfolder: report_FI_freq.py
    Input argument is a dictionary with a RePro indexes as a key and stimulus current as a value
    example: report_FI_freq.py -d"{'291': '-0.25nA', '293': '-0.25nA'}"
    """
    #   get the current working dir, split it into the dir and the path to it
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]
    #   read in some experimental data
    (_,_,filenames) = next(walk(ppjoin(wd, expfolder)))
    if "info.dat" in filenames:
        exp_info = dict(read_info(wd, expfolder)[0])
        print(exp_info)
    else:
        exp_info ={"Cell":{"Location":"UNLABELED"}}
    ReProList = command_interpreter(sys.argv[1:])
    print(ReProList["di"])
    fi_raster_return (ReProList["di"], exp_info)


