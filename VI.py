__author__ = 'jpresern'


import sys
sys.path.append("../pyrelacs")
# from IPython import embed
import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from itertools import groupby
from utils import seek_unique
from collections import OrderedDict

def vi_traces_sort(metas, datas):

    """
    Parses traces single VI curve RePro iteration and return them in a dictionary, where keys correspond to the
    stimuli used.
    :param metas : metas belonging to the selected RePro iteration
    :param datas : data belonging to the selected RePro iteration
    :return traceDict : spikes sorted together by the stimulus. A dictionary, of course.
    :return stim_amps : list of stimulus amplitudes as float. All values are assumed to be in nA.
    """

    StimList=[]
    for j in range(0, len(metas)):
        # extract the stimulus values
            StimList.append(metas[j]["I"])

    # StimIxCount = [len(list(group)) for key, group in groupby(StimList)]
    #   extracts the unique values only
    StimList = seek_unique(StimList)
    stim_amps = [float(x.strip("nA")) for x in StimList]
    traceDict = OrderedDict()
    # spikes_str_stimuli = {}
    for k, stimulus in enumerate(StimList):
        spikesIter = []
        for stimRep in range(len(metas)):
            if metas[stimRep]["I"] == stimulus:
                spikesIter.append(datas[stimRep])
        # print(spikesIter)
        traceDict[stimulus] = spikesIter

    return (traceDict, stim_amps)

def vi_curve(path_to_folder, subfolder, info):


    """"
    Plots the VI curve into the graphic file (*.png, *.svg)
    :param path_to_folder: location of the experimental folder
    :param subfolder: name of the experimental folder

    :return fnames: dict with VI figure filenames
    """
    #   define the input data
    filename1 = ppjoin(path_to_folder, subfolder, "vicurve-data.dat")
    filename2 = ppjoin(path_to_folder, subfolder, "vicurve-trace.dat")

    #   get the file, extract the three data subunits
    relacs_file1 = load(filename1)
    relacs_file2 = load(filename2)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        metas1, _, datas1 = relacs_file1.selectall()    #   get data from vicurve-data
        metas2, _, datas2 = relacs_file2.selectall()    #   get data from vicurve-trace

    except:
        return None

    ReProIxList = []
    for i in range(0,len(metas2)):
        comments2 = metas2[i]
        # extract the unique values only one time
        ReProIxList.append(comments2["ReProIndex"])

    # count occurences of same ReProIx
    ReProIx = [len(list(group)) for key, group in groupby(ReProIxList)]
    ReProIx.insert(0,0) # insert zero for easier job
    # print(ReProIx)

    #   count the iterations
    vi = np.array(datas1)
    VIiter = len(datas1)
    print(VIiter)

    #   setup plot dimension, 8" width, each RePro repetition gets 3" height
    # TODO: test figure output with a higher dpi settings (220 for retina display, 300 for print)
    # FHandles = plt.figure(figsize=(6, 3*VIiter))#, dpi=220)

    #   empty list for filename output
    fnames = []

    #   color scheme
    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "Dodgerblue", "Violet", "Silver", "Black", "Green", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "Dodgerblue", "Violet", "Silver", "Black", "Green" ]

    #   prepares the counter
    counter = np.cumsum(ReProIx)
    #   loop over VI RePro iterations
    for i in range(0, len(counter)-1):
        FHandles = plt.figure(figsize=(12, 3))#, dpi=220)
        print(i)
        comments_vi = metas1[i]

        #   counter increase
        #  counter = counter + ReProIx[i]
        #   convert the loaded data structures into an appropriate format

        #   sorts traces
        traceDict, stim_amps = vi_traces_sort(metas2[counter[i]:counter[i+1]], datas2[counter[i]:counter[i+1]])

        #   defines VI rates axis, number and position
        axarr = FHandles.add_subplot(1, 3, 1)

        #   draws plot
        # TODO: fix the errorbar color
        axarr.errorbar(vi[i][:, 0], vi[i][:, 7], vi[i][:, 8], color='red', fmt='-o')        #   peak
        axarr.errorbar(vi[i][:, 0], vi[i][:, 3], vi[i][:, 4], color='cyan', fmt='-x')       #   rest
        axarr.errorbar(vi[i][:, 0], vi[i][:, 5], vi[i][:, 6], color='yellow', fmt='-+')     #   steady state
        axarr.errorbar(vi[i][:, 0], vi[i][:, 10], vi[i][:, 11], color='green', fmt='-')     #   onset peak

        #   plot comments
        axarr.text(0.05, 0.95, " ".join(['Species:', info["Subject"]["Species"]]), transform=axarr.transAxes, fontsize = 10)
        axarr.text(0.05, 0.90, " ".join(['ELL Segment:', info["Cell"]["Location"]]), transform=axarr.transAxes, fontsize = 10)
        axarr.text(0.05, 0.85, " ".join(['RePro Ix:', str(comments_vi["ReProIndex"])]), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.80, " ".join(['RePro time:', str(comments_vi["ReProTime"])]), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.75, " ".join(['V_rest:', str(np.round(np.average(vi[i][:, 3]))), 'mV']), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.70, " ".join(['DC current:', str(vi[i][0, 1]), 'nA']), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.65, " ".join(['Amp. mode:', comments_vi["Status"]["AmplifierMode"]]), transform=axarr.transAxes, fontsize=10)
        if "SyncPulse" in comments_vi["Status"].keys():
            axarr.text(0.05, 0.60, " ".join(['DynClamp:', comments_vi["Status"]["SyncPulse"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.40, " ".join(['g:', comments_vi["Status"]["g"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.35, " ".join(['E:', comments_vi["Status"]["E"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.30, " ".join(['gVgate:', comments_vi["Status"]["gvgate"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.25, " ".join(['EVgate:', comments_vi["Status"]["Evgate"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.20, " ".join(['VgateTau:', comments_vi["Status"]["vgatetau"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.15, " ".join(['VgateMid:', comments_vi["Status"]["vgatevmid"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.10, " ".join(['VgateSlope:', comments_vi["Status"]["vgateslope"]]), transform=axarr.transAxes, fontsize=10)

        #   add traces subplot
        axarrM = FHandles.add_subplot(1,3,2)

        #   add resistance subplot
        axarrR = FHandles.add_subplot(1,3,3)

        #   compute the resistance
        resist = abs((vi[i][:, 5]-vi[i][:, 3])/vi[i][:,0])

        #   draws traces subplot and resistance subplot
        for j, k in enumerate(traceDict):
            axarrM.plot(traceDict[k][0][:,0],traceDict[k][0][:,1], color = cmap[j])
            axarrM.fill_between(traceDict[k][0][:,0], traceDict[k][0][:,1]-traceDict[k][0][:,2], traceDict[k][0][:,1]+traceDict[k][0][:,2], color = cmap[j], alpha=0.2)
            #   plot the resistance
            axarrR.plot(vi[i][j,0], resist[j], color = cmap[j], marker='o')

        #   writes titles, only over the top figures
        # if i == 0:
        axarr.set_title(subfolder + ": " + "VI curve", fontsize=12)

        #   turns off x axis labels for all except the bottom plot
        # if i < VIiter-1:
        # axarr.set_xticklabels([], visible=False)

        #   writes x labels beside the bottom most plot
        # if i == VIiter-1:
        axarr.set_xlabel('Current [nA]')
        axarrM.set_xlabel("time [ms]")
        axarrR.set_xlabel('Current [nA]')

        #   writes y labels
        axarr.set_ylabel('Voltage [mV]')
        axarrR.set_ylabel("Resistance [MOhm]")
        axarrM.set_yticklabels([], visible=False)

        #   define file name
        filename = "".join([str(comments_vi["ReProIndex"]),'_', 'VI_curve'])
        fnames.append(filename)

        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)
        plt.close()

    return fnames