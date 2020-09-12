__author__ = 'jpresern'

import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from itertools import groupby
from utils import seek_unique
from collections import OrderedDict
from FI import fi_spikes_sort, fi_instant_freq, fi_inst_freq_continuous

def sp_raster_plot(counter, path_to_folder, subfolder, spikeDict, freq_continuous, freq_continuous2,\
                timeDict, timeLine, metas, info, fileN):

    """
    sp_raster_plot plots "raster plots" for each stimulus and each stimulus iteration, for each stimulus iteration
    :param counter: which iteration of the RePro is being processes (important for file name)
    :param path_to_folder: location of the experimental folders
    :param subfolder: name of the experimental folder
    :param spikeDict: dictionary containing spikes for each stimulus iteration. key is the stimulus value
    :param metas: meta info from relacs files
    :param freq_continuous: instantaneous frequency computed for each time point in timeLine for all stimulus iterations
    :param freq_continuous2: instantaneous frequency computed fore ech time point in timeLine for all stimulus
            iterations which respond with spiking during stimulus presentation
    """
    #   defines relative sizes of individual raster plots
    rug_bot = 0.05
    rug_width  = 0.60
    rug_height = 0.67/len(spikeDict)
    rug_step = 1/len(spikeDict)

    #   defines relative size of individual return plots
    return_bot = 0.70
    return_width = 0.1
    return_height = 0.67/len(spikeDict)
    return_step = 1/len(spikeDict)

    #   defines relative size of individual histograms
    hist_bot =  0.88
    hist_width = 0.1
    hist_height = 0.67/len(spikeDict)
    hist_step = 1/len(spikeDict)

    #   define the colormap
    cmapBlue = ["DarkTurquoise","DarkCyan", "CadetBlue", "DeepSkyBlue", "CornFlowerBlue", "DodgerBlue", "LightSkyBlue", "LightSteelBlue"]
    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]
    #   defines absolute size of the rug figure (considering the number of different stimuli
    FHandles = plt.figure(figsize=(15, 2*len(spikeDict)))
    for i, k in enumerate(spikeDict):

        axarrL = FHandles.add_axes ([rug_bot,(i)*rug_step+(0.25/len(spikeDict)), rug_width, rug_height])
        axarrL2 = axarrL.twinx()

        axarrL.plot(timeLine, freq_continuous[k], color="g")
        axarrL.plot(timeLine, freq_continuous2[k], color="b")

        # axarrL2.plot(fitTimeLine, instFreqFit[k], color="r")

        for j in range(len(spikeDict[k])):
            axarrL2.plot(spikeDict[k][j], np.zeros(len(spikeDict[k][j]))+j, '|', color='red', ms= 12)

        axarrL.text(0.05, 0.90, " ".join(['Species:', info["Subject"]["Species"]]), transform=axarrL.transAxes, fontsize = 10)
        axarrL.text(0.05, 0.80, " ".join(['Stimulus:', k]), transform=axarrL.transAxes, fontsize=10)
        axarrL.text(0.05, 0.70, " ".join(['DC current:', metas["Status"]["Current-1"]]), transform=axarrL.transAxes, fontsize=10)
        # axarrL.text(0.05, 0.75, " ".join(['No. Repro iter:', metas[i]["trials"]]), transform=axarrL.transAxes, fontsize=10)

        if "SyncPulse" in metas["Status"].keys():
            axarrL.text(0.40, 0.90, " ".join(['DynClamp:', metas["Status"]["SyncPulse"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.80, " ".join(['g:', metas["Status"]["g"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.70, " ".join(['E:', metas["Status"]["E"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.90, " ".join(['gVgate:', metas["Status"]["gvgate"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.80, " ".join(['EVgate:', metas["Status"]["Evgate"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.70, " ".join(['VgateTau:', metas["Status"]["vgatetau"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.60, " ".join(['VgateMid:', metas["Status"]["vgatevmid"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.50, " ".join(['VgateSlope:', metas["Status"]["vgateslope"]]), transform=axarrL.transAxes, fontsize=10)

        axarrL.set_xlim(0, float(metas["Settings"]["General"]["duration"].split('s')[0])*1000)
        # axarrL.set_xticklabels([0,5,10,15,20,25,30])
        axarrL2.set_yticklabels([], visible=False)
        axarrL2.set_ylim(0,len(spikeDict[k]))

        #   add y labels
        # axarrL2.set_ylabel("Stim Iter")
        axarrL.set_ylabel("Freq [Hz]", color = "b")
        axarrL.set_xlabel("Time [ms]")
        #   paint the axis red
        for tl in axarrL.get_yticklabels():
            tl.set_color("b")

        #   drawing the second plot, containing time intervals
        axarrM = FHandles.add_axes([return_bot,(i)*return_step+(0.25/len(spikeDict)), return_width, return_height])

        allIntervalsZero = np.array([0])
        for l in range(len(timeDict[k])):
            allIntervals = np.append(allIntervalsZero, timeDict[k][l])

            for m in range(len(allIntervals)-1):
               axarrM.plot(allIntervals[m],allIntervals[m+1], "o", color= cmap[l], ms=5)

        # axarrR.set_xlim(0, np.percentile(allIntervals,97))
        # axarrR.set_ylim(0, np.percentile(allIntervals,97))

        axarrM.set_xlim(0, 35)
        axarrM.set_ylim(0, 35)

        axarrM.set_ylabel("ISI [ms]")
        axarrM.set_xlabel("ISI [ms]")

        #   drawing the third plot, containing time interval distribution
        axarrR = FHandles.add_axes([hist_bot,(i)*hist_step+(0.25/len(spikeDict)), hist_width, hist_height])
        axarrR.set_ylabel("fraction")
        axarrR.set_xlabel("ISI [ms]")

        #   put all insterspike intervals into a single numpy array
        allIntervals = []
        for l in range(len(timeDict[k])):
            allIntervals.extend(timeDict[k][l])
        allIntervals = np.array(allIntervals)

        #   histogram calculation for the ISI distribution
        if allIntervals.shape[0] > 1:
            allIntervals = allIntervals[allIntervals<np.percentile(allIntervals,95)]
            if allIntervals.shape[0] > 1:
                axarrR.hist (allIntervals, bins = 29, normed = True)

    #   define file name
    filename = "".join([str(fileN),'_', 'spont_act'])

    FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
    FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)

    # plt.show()
    plt.close()

    return filename


def spont_act(path_to_folder, subfolder, info):

    """
    main calls various sub-functions to exctract the data and plot the raster, return map and interval histogram
    :param path_to_folder: folder containing experiments
    :param subfolder: folder containing cell experimental data
    :param info: experiment metadata
    :return fnames
    """
    #   define the input data
    filename = ppjoin(path_to_folder, subfolder, "saveevents-Spikes-1.dat")

    #   get the file, extract the three data subunits
    relacs_file = load(filename)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        metas, _, datas = relacs_file.selectall()

    except:
        return None

    #   define list for rug filenames
    fnames_spont = []

    #   for each iteration of the same RePro
    for i in range(0, len(metas)):

        #   print processed RePro iteration
        print("Spont activity ReProIx", i)

        #   sorts spikes
        aa = [metas[i]]
        #   conversion into miliseconds
        spikeDict, stim_amps = fi_spikes_sort(aa, datas[i][:,0]*1000)

        #   computes instantaneous spike frequencies based on spikes
        freqDict, timeDict = fi_instant_freq(spikeDict) #, metas[counter[i]:counter[i+1]])
        freq_continuous, freq_continuous2, timeLine, freqSD_continuous = fi_inst_freq_continuous(spikeDict, freqDict, metas)

        #   plots
        fnames = sp_raster_plot(i, path_to_folder, subfolder, spikeDict, freq_continuous, freq_continuous2, timeDict,\
                    timeLine, metas[i], info, fileN = i)

        #   appends the list of figure file names
        fnames_spont.append(fnames)

    return fnames_spont
