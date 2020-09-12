__author__ = 'jpresern'

import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
# from os.path import join
from posixpath import join as ppjoin
from itertools import groupby
from utils import seek_unique
from collections import OrderedDict
from mathEngine import exp_decay2

# TODO: extract first spike latency, plot it against current stimulus
# TODO: fit exponential decay over the instantaneous frequencies
# TODO: plot time constants against stimulus.

def fi_spikes_sort(metas, datas):

    """
    Parses spikes from single FI curve RePro iteration and return them in a dictionary, where keys correspond to the
    stimuli used.
    :param metas: metas belonging to the selected RePro iteration
    :param datas: data belonging to the selected RePro iteration
    :return spikes_str_stimuli: spikes sorted together by the stimulus. A dictionary, of course.
    :return stim_amps: list of stimulus amplitudes as float. All values are assumed to be in nA.
    """

    StimList=[]
    for j in range(0, len(metas)):
        # extract the stimulus values
            try:
                StimList.append(metas[j]["I"])

            except:
                StimList.append("-0.0nA")

    # StimIxCount = [len(list(group)) for key, group in groupby(StimList)]
    #   extracts the unique values only
    StimList = seek_unique(StimList)
    stim_amps = [float(x.strip("nA")) for x in StimList]
    spikes_str_stimuli = OrderedDict()
    # spikes_str_stimuli = {}
    for k, stimulus in enumerate(StimList):
        spikesIter = []
        for stimRep in range(len(metas)):
            #   checks if the metas contain abbreviated (corrupted?) header (saveevents-Spikes-1.dat)
            if "events" in metas[stimRep]:
                if metas[stimRep]["events"] == "Spikes-1":
                    spikesIter.append(datas)
            elif "I" in metas[stimRep]:
                if metas[stimRep]["I"] == stimulus:
                    spikesIter.append(datas[stimRep])

        # print(spikesIter)
        spikes_str_stimuli[stimulus] = spikesIter

    return (spikes_str_stimuli, stim_amps)

def fi_instant_freq(spikeDict):

    """
    Parses dictionary, containing spike (times) and calculates the instantaneous frequency
    :param spikeDict : a dictionary, containing stimuli as keys and list of lists for each entry. Spikes are
                        written as nd.array in the bottom list
    :param metas : meta containing duration of stimulus, delay and pause.

    :returns freqDict : dictionary containing instant spike frequencies
    :returns timeDict : dictionary containing time differences between spikes
    """
    freqDict = OrderedDict()
    timeDict = OrderedDict()
    for k, v in spikeDict.items():
        freqIter = []
        timeIter = []
        for i in range(len(v)):
            timeDiff = np.diff(v[i])
            inst_freq = (1/timeDiff)*1000
            freqIter.append(inst_freq)
            timeIter.append(timeDiff)
        freqDict[k] = freqIter
        timeDict[k] = timeIter
    return freqDict, timeDict

def fi_inst_freq_continuous(spikeDict, freqDict, metas):

    """
    fi_inst_freq_continuous processes instantaneous spike frequency according to event. Function constructs
    a dictionary of lists of np.arrays, which are actually continuous vectors of spike frequencies.
    It assumes that freq remains unchanged between two consecutive spikes.
    : param spikeDict : ordered dictionary of spikes, arranged by stimuli and iterations of same stimuli
    : param freqDict : ordered dictionary of frequencies between spikes, n-1 values, arranged by stimuli and iterations of same stimuli
    : param metas : meta info containing, from which duration of single stimulus iteration is extracted.
    : return freq_continuous : an ordered dictionary of lists of np.array, containing averaged frequency values
    : return freq_continuous2 : an ordered dictionary of lists of np.array, containing averaged frequency values of inst. freq. during stimulus only
    : return timeLine : a np.array of time points
    """

    #   definining time line: if dt = 1 ms, number of points equals duration in ms
    dt = 1 #    1/samplerate in miliseconds
    #   checks if the metas contain abbreviated (corrupted?) header (saveevents-Spikes-1.dat)
    if "events" in metas[0]:
        if metas[0]["events"] == "Spikes-1":
            duration = float(metas[0]["Settings"]["General"]["duration"].split('s')[0])*1000
            timeLine = np.linspace(0, duration, (duration+1)/dt)
    else:
        duration = float(metas[0]["Settings"]["Timing"]["delay"].split('m')[0]) \
                        + float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0])+300

        #   defining timeLine: it is duration of the delay and stimulus extended for 300 ms. An additional ms is added
        #   to cover for 0.
        timeLine = np.linspace(-float(metas[0]["Settings"]["Timing"]["delay"].split('m')[0]),\
                                float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0])+300, (duration+1)/dt)

    #   constructing OrderedDictionary to contain averaged instantaneous frequencies
    freq_continuous = OrderedDict()
    freqSD_continuous = OrderedDict()
    #   constructing OrderedDictionary to contain averaged instantaneous frequencies that respond with spikes during stimuli
    freq_continuous2 = OrderedDict()
    #   looping over spikeDict/freqDict
    for i, k in enumerate(spikeDict):
        #   constructing the frequency vector corresponding to the timeLine
        freq = np.zeros(shape = (len(spikeDict[k]),len(timeLine)))
        freqSD = np.zeros(shape = (len(spikeDict[k]),len(timeLine)))
        # print(k)
        for j in range(len(freqDict[k])):
            #   counter remembers which of the spike trains were empty
            counter = []
            #   counter2 remembers which of the spike trains were empty during the stimulus time
            counter2 = []
            #   cheks if there is an empty array
            #   TODO: add option for the corrupted (abbreviated) header
            if (freqDict[k][j].shape != (0,)) & (spikeDict[k][j].any() != -0.):
                # freq = np.row_stack(np.zeros(shape = (1,len(timeLine))))
                #   replaces zeros by interspike intervals
                for l,spike in enumerate(freqDict[k][j]):
                    freq[j, (timeLine>=spikeDict[k][j][l]) & (timeLine<spikeDict[k][j][l+1])] = freqDict[k][j][l]

                # if (0<=spikeDict[k][j]) & (spikeDict[k][j]<float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0])) == []:
                    # counter2.append(j)
                    # print("v stimulus ni nic", k, j)

            #   count how many lines were empty
            # else:
                # counter.append(j)
                # print("nema nista", k, j)

        ## remove empty lines before averaging
        # freq2 = np.delete(freq, counter2, 0)
        # freq = np.delete(freq, counter, 0)
        #   calculates mean over columns
        if freq.shape[0] >= 1:
            freqSD = np.std(freq,0)
            freq = np.mean(freq, 0)

        elif freq.shape[0] == 0:
            freq = np.zeros(shape = (len(timeLine)))

        # if freq2.shape[0] >= 1:
        #     freq2 = np.mean(freq2, 0)
        # elif freq2.shape[0] == 0:
        #     freq2 = np.zeros(shape = (len(timeLine)))

        #   insert computed frequency into the OrderedDict
        freq_continuous[k] = freq
        freqSD_continuous[k] = freqSD
        # freq_continuous2[k] = freq2
        freq_continuous2[k] = freq
    return freq_continuous, freq_continuous2, timeLine, freqSD_continuous

def fi_tau_fit (freq_continuous, metas):

    """
    A function that fits exponential decay over the averaged instantaneous frequency.
    : param fi_isnt_freq_continuous : an ordered dictionary, with list of np.arrays for each trial
    : return taus : an ordered dictionary with taus for each stimulus
    : return instFreqFit : an ordered dictionary with computed values for the fit
    : return fitTimeLine : a numpy.array containing the corresponding time values used in the fit.
    """
    #   definining time line: if dt = 1 ms, number of points equals duration in ms
    dt = 1 #    1/samplerate in miliseconds
    duration = float(metas[0]["Settings"]["Timing"]["delay"].split('m')[0]) \
                        + float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0])+300

    #   defining stim duration
    stimDur = float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0])

    #   defining timeLine: it is duration of the delay and stimulus extended for 300 ms. An additional ms is added
    #   to cover for 0.
    timeLine = np.linspace(-float(metas[0]["Settings"]["Timing"]["delay"].split('m')[0]),\
                            stimDur+300, (duration+1)/dt)

    #   time line used in computing the fitted values
    fitTimeLine = timeLine[(timeLine > 0) & (timeLine <=stimDur)]

    taus = OrderedDict()
    instFreqFit = OrderedDict()
    for i, k in enumerate(freq_continuous):
        if freq_continuous[k][(timeLine > 0) & (timeLine <= stimDur)].all() != 0:

            #   fits the exponential decay over the voltage response to the current pulse
            try:
                parameters = exp_decay2(fitTimeLine,freq_continuous[k][(timeLine > 0) & (timeLine <=stimDur)],(max(freq_continuous[k][(timeLine > 0) & (timeLine <=stimDur)]), 9, np.mean(freq_continuous[k][(timeLine > 0.9*stimDur) & (timeLine <=stimDur)])))
            except:
                parameters = [0, 0, 0]
        else:
            parameters = [0, 0, 0]


        freqFit = parameters[0]*np.exp(-parameters[1]*fitTimeLine)\
                     +parameters[2]

        print(k, parameters)
        taus[k] = parameters
        instFreqFit[k] = freqFit

    return taus, instFreqFit, fitTimeLine

def dyn_fi_freq(RePro, wd, expfolder):

    """
    dyn_fi_freq used FI.fi_spikes_sort, FI.fi_instant_freq and FI.fi_inst_freq_continuous to extract the steady state and peak frequencies
    :param RePro: list of RePros to be processes
    :param wd: working directory
    :param expfolder: folder containing actual experiments
    :return ss: ordered dictionary containing a list under ReProIndex. This list contains 1) nd array with desired frequencies paired with currents, gvgate and vgatetau and g values.
    :return peak: ordered dictionary containing a list under ReProIndex. This list contains 1) nd array with desired frequencies paired with currents, gvgate and vgatetau and g values.
    :return rest: ordered dictionary containing a list under ReProIndex. This list contains 1) nd array with desired frequencies paired with currents, gvgate and vgatetau and g values.
    :return metas: metas
    """

    #   define sample rate (instantaneous frequenices are computed with the 1 ms resolution)
    dt = 1

    #   define top dictionaries
    ss = OrderedDict()
    peak = OrderedDict()
    rest = OrderedDict()

    for ix in RePro:

        #   define the input data
        filename = ppjoin(wd, expfolder, "ficurve-spikes.dat")

        #   get the file, extract the three data subunits
        relacs_file = load(filename)

        #   if relacs file is empty or too short due to aborted RePro presentation
        #   "try:" provides normal termination of the loop instead of error
        try:
            metas, _, datas = relacs_file.select({"ReProIndex": ix})
        except:
            return None

        #   extracts RePro indexes
        ReProIxList = []
        for i in range(0,len(metas)):
            # comments2 = metas[i]
            # extract the unique values only one time
            ReProIxList.append(metas[i]["ReProIndex"])

        # count occurrences of same ReProIx
        ReProIx = [len(list(group)) for key, group in groupby(ReProIxList)]

        ReProIx.insert(0, 0) # insert zero for easier job in indexing
        counter = np.cumsum(ReProIx)

        #   define list for rug filenames
        fnames_rug = []

        #   print processed RePro iteration
        print("FI Raster & Return ReProIx", ReProIxList[counter[0]])#metas[ReProIx[i]]["ReProIndex"])

        #   sorts spikes
        spikeDict, stim_amps = fi_spikes_sort(metas[counter[0]:counter[1]], datas[counter[0]:counter[1]])
        #   computes instantaneous spike frequencies based on spikes
        freqDict, timeDict = fi_instant_freq(spikeDict) #, metas[counter[i]:counter[i+1]])
        freq_continuous, freq_continuous2, timeLine, freqSD_continuous = fi_inst_freq_continuous(spikeDict, freqDict, metas)

        #   define empty nd arrays to contain steady states and peak freqs alongside with injected currents
        steady_state = np.ndarray(shape=(len(freq_continuous),3), dtype=float)
        peak_freq = np.ndarray(shape=(len(freq_continuous),3), dtype=float)
        rest_freq = np.ndarray(shape=(len(freq_continuous),3), dtype=float)

        #   extract duration of initial interval without stimulus and duration of stimulus
        durDelay = int(metas[0]["Settings"]["Timing"]["delay"].strip("ms"))
        durStimulus = int(metas[0]["Settings"]["Timing"]["duration"].strip("ms"))

        counter2 = 0

        for k, v in freq_continuous.items():
        #   TODO: automatically obtain test pulse length
            steady_state[counter2,0] = float(k.strip("nA"))
            steady_state[counter2,1] = np.average(freq_continuous[k][durStimulus+durDelay-150:durStimulus+durDelay])
            steady_state[counter2,2] = np.average(freqSD_continuous[k][durStimulus+durDelay-150:durStimulus+durDelay])

            peak_freq[counter2,0] = float(k.strip("nA"))
            # peak_freq[counter2,1] = np.max(freq_continuous[k][durDelay:durDelay+100])   #   computes peak freq in certain time window
            peak_freq[counter2,1] = np.max(freq_continuous[k][durDelay:durDelay+durStimulus])   #   takes out maximum frequency event
            peak_freq[counter2,2] = freqSD_continuous[k][durDelay - 1 + np.argmax(freq_continuous[k][durDelay:durDelay+100])]

            rest_freq[counter2,0] = float(k.strip("nA"))
            rest_freq[counter2,1] = np.max(freq_continuous[k][0:durDelay])
            rest_freq[counter2,2] = np.average(freqSD_continuous[k][0:durDelay])

            counter2 = counter2+1

        ss[ix] = [steady_state, metas[0]["Status"]["gvgate"], metas[0]["Status"]["vgatetau"], metas[0]["Status"]["g"]]
        peak[ix] = [peak_freq, metas[0]["Status"]["gvgate"], metas[0]["Status"]["vgatetau"], metas[0]["Status"]["g"]]
        rest[ix] = [rest_freq, metas[0]["Status"]["gvgate"], metas[0]["Status"]["vgatetau"], metas[0]["Status"]["g"]]

    return ss, peak, metas, rest

def fi_rug_plot(counter, path_to_folder, subfolder, spikeDict, freq_continuous, freq_continuous2,\
                timeDict, timeLine, metas, info, fileN = 'FI_rug'):

    """
    fi_rug_plot plots "raster plots" for each stimulus and each stimulus iteration, for each stimulus iteration
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
    rug_width  = 0.48
    rug_height = 0.67/len(spikeDict)
    rug_step = 1/len(spikeDict)

    #   defines relative size of individual return plots
    return_bot = 0.63
    return_width = 0.13
    return_height = 0.67/len(spikeDict)
    return_step = 1/len(spikeDict)

    #   defines relative size of individual histograms
    hist_bot =  0.83
    hist_width = 0.13
    hist_height = 0.67/len(spikeDict)
    hist_step = 1/len(spikeDict)


    #   define the colormap
    cmapBlue = ["DarkTurquoise","DarkCyan", "CadetBlue", "DeepSkyBlue", "CornFlowerBlue", "DodgerBlue", "LightSkyBlue", "LightSteelBlue"]
    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]
    #   defines absolute size of the rug figure (considering the number of different stimuli
    FHandles = plt.figure(figsize=(14, 2*len(spikeDict)))
    for i, k in enumerate(spikeDict):


        axarrL = FHandles.add_axes ([rug_bot,(i)*rug_step+(0.25/len(spikeDict)), rug_width, rug_height])
        axarrL2 = axarrL.twinx()

        axarrL2.plot(timeLine, freq_continuous[k], color="g")
        axarrL2.plot(timeLine, freq_continuous2[k], color="b")

        # axarrL2.plot(fitTimeLine, instFreqFit[k], color="r")

        for j in range(len(spikeDict[k])):
            axarrL.plot(spikeDict[k][j], np.ones(len(spikeDict[k][j]))+j, '|', color=cmap[j], ms= 12)

        axarrL.text(0.05, 0.90, " ".join(['Species:', info["Subject"]["Species"]]), transform=axarrL.transAxes, fontsize = 10)
        axarrL.text(0.05, 0.80, " ".join(['Stimulus:', k]), transform=axarrL.transAxes, fontsize=10)
        axarrL.text(0.05, 0.70, " ".join(['DC current:', metas[0]["Status"]["Current-1"]]), transform=axarrL.transAxes, fontsize=10)
        # axarrL.text(0.05, 0.75, " ".join(['No. Repro iter:', metas[i]["trials"]]), transform=axarrL.transAxes, fontsize=10)

        if "SyncPulse" in metas[0]["Status"].keys():
            axarrL.text(0.70, 0.70, " ".join(['DynClamp:', metas[0]["Status"]["SyncPulse"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.60, " ".join(['g:', metas[0]["Status"]["g"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.50, " ".join(['E:', metas[0]["Status"]["E"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.40, " ".join(['gVgate:', metas[0]["Status"]["gvgate"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.30, " ".join(['EVgate:', metas[0]["Status"]["Evgate"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.20, " ".join(['VgateTau:', metas[0]["Status"]["vgatetau"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.10, " ".join(['VgateMid:', metas[0]["Status"]["vgatevmid"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.70, 0.0, " ".join(['VgateSlope:', metas[0]["Status"]["vgateslope"]]), transform=axarrL.transAxes, fontsize=10)


        axarrL.set_xlim(-float(metas[0]["Settings"]["Timing"]["delay"].split('m')[0]),\
                        float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0])+300)
        axarrL.set_ylim(0,len(spikeDict[k])+1)

        axarrL.axvline(0., color='y')
        axarrL.axvline(float(metas[0]["Settings"]["Timing"]["duration"].split('m')[0]), color='y')

        #   add y labels
        axarrL.set_ylabel("Stim Iter")
        axarrL2.set_ylabel("Freq [Hz]", color = "b")
        #   paint the axis blue
        for tl in axarrL2.get_yticklabels():
            tl.set_color("b")

        #   drawing the second plot, containing time intervals
        axarrM = FHandles.add_axes([return_bot,(i)*return_step+(0.25/len(spikeDict)), return_width, return_height])

        allIntervalsZero = np.array([0])
        for l in range(len(timeDict[k])):
            allIntervals = np.append(allIntervalsZero, timeDict[k][l])

            for m in range(len(allIntervals)-1):
               axarrM.plot(allIntervals[m],allIntervals[m+1]-allIntervals[m], "o", color= cmap[l], ms=5)

        axarrM.set_xlim(0, 14)
        axarrM.set_ylim(-7, 7)

        #   drawing the third plot, containing time interval distribution
        axarrR = FHandles.add_axes([hist_bot,(i)*hist_step+(0.25/len(spikeDict)), hist_width, hist_height])

        #   put all insterspike intervals into a single numpy array
        allIntervals = []
        for znj in range(len(timeDict[k])):
            IntervalsN = np.append(allIntervalsZero, timeDict[k][znj])
            # for ml in range((IntervalsN).shape[0]-1):
            #     allIntervals.append((IntervalsN[ml+1]-IntervalsN[ml]))

            IntervalsN1 = np.append(timeDict[k][znj], allIntervalsZero)
            allIntervals.extend((IntervalsN1-IntervalsN)[1:-2])

        allIntervals = np.array(allIntervals)

        axarrM.set_xlabel('n ISI [ms]')
        axarrM.set_ylabel('n+1 ISI - n ISI [ms]')

        # allIntervals = []
        # for l in range(len(timeDict[k])):
        #     allIntervals.extend(timeDict[k][l])
        # allIntervals = np.array(allIntervals)

        #   histogram calculation for the ISI distribution
        if allIntervals.shape[0] > 1:
            # allIntervals = allIntervals[allIntervals<np.percentile(allIntervals,95)]
            if allIntervals.shape[0] > 1:
                axarrR.hist (allIntervals, bins = 20, normed = False, orientation = 'horizontal')

        #   add labels
        axarrR.set_ylabel('n+1 ISI - n ISI [ms]')
        axarrR.set_xlabel('Count')
        axarrR.set_ylim(-7, 7)



    #   define file name
    filename = "".join([str(metas[0]["ReProIndex"]),'_', fileN])

    FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
    FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)

    # plt.show()
    plt.close()

    return filename

def fi_curve(path_to_folder, subfolder, info):

    """"
    Plots the FI curve as computed by Relacs into the graphic file (*.png, *.svg)
        Requires:
    :param path_to_folder: location of the experimental folder
    :param subfolder: name of the experimental folder
        Outputs:
            graphic files containing FI plots (rates and PSTH) as .png and .svg
    """
    #   define the input data
    filename1 = ppjoin(path_to_folder, subfolder, "ficurve-data.dat")
    filename2 = ppjoin(path_to_folder, subfolder, "ficurve-rates.dat")

    #   get the file, extract the three data subunits
    relacs_file1 = load(filename1)
    relacs_file2 = load(filename2)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        metas1, _, datas1 = relacs_file1.selectall()
        metas2, _, datas2 = relacs_file2.selectall()

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
    FIiter = len(datas1)
    print(FIiter)

     #   empty list for filename output
    fnames = []
    counter = 0
    for i in range(0,FIiter):

        #   print processed RePro iteration
        print("FI Curve ReProIx", metas1[i]["ReProIndex"])

        # TODO: test figure output with a higher dpi settings (220 for retina display, 300 for print)
        FHandles = plt.figure(figsize=(10, 4))#, dpi=220)
        #   counter increase
        counter = counter + ReProIx[i]
        #   convert the loaded data structures into an appropriate format
        rates1 = np.array(datas1[i])
        rates2 = np.array(datas2[counter: counter + ReProIx[i+1]])
        datas1[i] = None
        datas2[i] = None
        comments1 = metas1[i]

        #   FI frequencies - select data
        #   select the data columns

        Istep = rates1[:, 0]
        FRavg = rates1[:, 3]
        FRavgSD = rates1[:, 4]
        FRbase = rates1[:, 5]
        FRbaseSD = rates1[:, 6]
        FRpeak = rates1[:, 9]
        FRpeakSD = rates1[:, 10]
        FRsteady = rates1[:, 12]
        FRsteadySD = rates1[:, 13]

        #   defines FI rates axis, number and position
        axarrL = FHandles.add_subplot(1, 2, 1)

        #   plots FI curves into the axis
        axarrL.plot(Istep, FRbase, color='b')
        axarrL.fill_between(Istep, FRbase+FRbaseSD, FRbase-FRbaseSD, color='b', alpha=0.2)
        axarrL.plot(Istep, FRpeak, color='r')
        axarrL.fill_between(Istep, FRpeak+FRpeakSD, FRpeak-FRpeakSD, color='r', alpha=0.2)
        axarrL.plot(Istep, FRsteady, color='g')
        axarrL.fill_between(Istep, FRsteady+FRsteadySD, FRsteady-FRsteadySD, color='g', alpha=0.2)
        axarrL.plot(Istep, FRavg, color='m')
        axarrL.fill_between(Istep, FRavg+FRavgSD, FRavg-FRavgSD, color='m', alpha=0.2)

        #   plot comments
        axarrL.text(0.05, 0.95, " ".join(['Species:', info["Subject"]["Species"]]), transform=axarrL.transAxes, fontsize = 10)
        axarrL.text(0.05, 0.90, " ".join(['ELL Segment:', info["Cell"]["Location"]]), transform=axarrL.transAxes, fontsize = 10)
        axarrL.text(0.05, 0.85, " ".join(['RePro Ix:', str(comments1["ReProIndex"])]), transform=axarrL.transAxes, fontsize=10)
        axarrL.text(0.05, 0.80, " ".join(['RePro time:', str(comments1["ReProTime"])]), transform=axarrL.transAxes, fontsize=10)
        axarrL.text(0.05, 0.75, " ".join(['V_rest:', str(np.round(np.average(rates1[:, 7]))), 'mV']), transform=axarrL.transAxes, fontsize=10)
        axarrL.text(0.05, 0.70, " ".join(['DC current:', str(rates1[0, 1]), 'nA']), transform=axarrL.transAxes, fontsize=10)
        axarrL.text(0.05, 0.65, " ".join(['Amp. mode:', comments1["Status"]["AmplifierMode"]]), transform=axarrL.transAxes, fontsize=10)

        if "SyncPulse" in comments1["Status"].keys():
            axarrL.text(0.40, 0.60, " ".join(['SyncPulse:', comments1["Status"]["SyncPulse"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.40, " ".join(['g:', comments1["Status"]["g"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.35, " ".join(['E:', comments1["Status"]["E"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.30, " ".join(['gVgate:', comments1["Status"]["gvgate"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.25, " ".join(['EVgate:', comments1["Status"]["Evgate"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.20, " ".join(['VgateTau:', comments1["Status"]["vgatetau"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.15, " ".join(['VgateMid:', comments1["Status"]["vgatevmid"]]), transform=axarrL.transAxes, fontsize=10)
            axarrL.text(0.40, 0.10, " ".join(['VgateSlope:', comments1["Status"]["vgateslope"]]), transform=axarrL.transAxes, fontsize=10)


        #   defines FI PSTH axis, number and position
        axarrR = FHandles.add_subplot(1, 2, 2)

        #   plots FI PSTH
        for j in range(0, len(rates2)):
            axarrR.plot(rates2[j][:, 0]/1000, rates2[j][:, 1])#, line_color=colors[j])

        #   sets y label on the left plot
        axarrL.set_ylabel('Frequency [Hz]')
        axarrR.set_yticklabels([], visible=False)

        #   sets y lim on both plots
        axarrL.set_ylim(axarrR.get_ylim())

        #   writes titles, only over the top figures
        axarrL.set_title(subfolder + ": " + "FI rates", fontsize=12)
        axarrR.set_title(subfolder + ": " + "FI PSTH", fontsize=12)

        #   writes labels beside the bottom most plot
        axarrL.set_xlabel('Current [nA]')
        axarrR.set_xlabel('time [s]')

        #   define file name
        filename = "".join([str(metas1[i]["ReProIndex"]),'_', 'FI_curve'])
        fnames.append(filename)

        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)

        plt.close()

    return fnames

def fi_curve2(counter, path_to_folder, subfolder, timeLine, freq_continuous, freqSD_continuous, metas, info):

    """
    function fi_curve2 extracts spike frequencies at steady state, peak and rest and plots them
    :param path_to_folder: location of the experimental folder
    :param subfolder: experimental folder itself
    :param freq_continuous: time line of frequencies (a dict)
    :param freqSD_continuous: time line of SD of frequencies (a dict)
    :return fnames: list of file names of the plots
    """

 #   define empty nd arrays to contain steady states and peak freqs alongside with injected currents
    steady_state = np.ndarray(shape=(len(freq_continuous),3), dtype=float)
    peak_freq = np.ndarray(shape=(len(freq_continuous),3), dtype=float)
    rest_freq = np.ndarray(shape=(len(freq_continuous),3), dtype=float)

    #   extract duration of initial interval without stimulus and duration of stimulus
    durDelay = int(metas[0]["Settings"]["Timing"]["delay"].strip("ms"))
    durStimulus = int(metas[0]["Settings"]["Timing"]["duration"].strip("ms"))

    counter2 = 0

    #TODO: write experimental metas on the figures!

    #   define drawing area
    FHandles = plt.figure(figsize=(10, 4))
    axFI   = FHandles.add_axes([0.10, 0.15, 0.40, 0.77])
    axFreq = FHandles.add_axes([0.55, 0.15, 0.40, 0.77])

    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]

    for k, v in freq_continuous.items():
    #   TODO: automatically obtain test pulse length
        steady_state[counter2,0] = float(k.strip("nA"))
        steady_state[counter2,1] = np.average(freq_continuous[k][durStimulus+durDelay-150:durStimulus+durDelay])
        steady_state[counter2,2] = np.average(freqSD_continuous[k][durStimulus+durDelay-150:durStimulus+durDelay])

        peak_freq[counter2,0] = float(k.strip("nA"))
        peak_freq[counter2,1] = np.max(freq_continuous[k][durDelay:durDelay+100])
        peak_freq[counter2,2] = freqSD_continuous[k][durDelay - 1 + np.argmax(freq_continuous[k][durDelay:durDelay+100])]

        rest_freq[counter2,0] = float(k.strip("nA"))
        rest_freq[counter2,1] = np.max(freq_continuous[k][0:durDelay])
        rest_freq[counter2,2] = np.average(freqSD_continuous[k][0:durDelay])

        axFreq.plot(timeLine, freq_continuous[k], color=cmap[counter2])

        counter2 = counter2+1

    axFI.errorbar(steady_state[:,0], steady_state[:,1], steady_state[:,2], color="Gold", marker="o")
    axFI.errorbar(peak_freq[:,0], peak_freq[:,1], peak_freq[:,2], color="Red", marker="o")
    axFI.errorbar(rest_freq[:,0], rest_freq[:,1], rest_freq[:,2], color="DodgerBlue", marker="o")

    #   plot comments
    axFI.text(0.05, 0.95, " ".join(['Species:', info["Subject"]["Species"]]), transform=axFI.transAxes, fontsize = 10)
    axFI.text(0.05, 0.90, " ".join(['ELL Segment:', info["Cell"]["Location"]]), transform=axFI.transAxes, fontsize = 10)
    axFI.text(0.05, 0.85, " ".join(['RePro Ix:', str(metas[counter]["ReProIndex"])]), transform=axFI.transAxes, fontsize=10)
    axFI.text(0.05, 0.80, " ".join(['RePro time:', str(metas[counter]["ReProTime"])]), transform=axFI.transAxes, fontsize=10)
    # axFI.text(0.05, 0.75, " ".join(['V_rest:', str(np.round(np.average(rates1[:, 7]))), 'mV']), transform=axFI.transAxes, fontsize=10)
    axFI.text(0.05, 0.70, " ".join(['DC current:', str(metas[counter]["Status"]["Current-1"])]), transform=axFI.transAxes, fontsize=10)
    axFI.text(0.05, 0.65, " ".join(['Amp. mode:', metas[counter]["Status"]["AmplifierMode"]]), transform=axFI.transAxes, fontsize=10)


    if "SyncPulse" in metas[counter]["Status"].keys():

        axFI.text(0.05, 0.50, " ".join(['SyncPulse:', metas[counter]["Status"]["SyncPulse"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.45, " ".join(['g:', metas[counter]["Status"]["g"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.40, " ".join(['E:', metas[counter]["Status"]["E"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.35, " ".join(['gVgate:', metas[counter]["Status"]["gvgate"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.30, " ".join(['EVgate:', metas[counter]["Status"]["Evgate"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.25, " ".join(['VgateTau:', metas[counter]["Status"]["vgatetau"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.20, " ".join(['VgateMid:', metas[counter]["Status"]["vgatevmid"]]), transform=axFI.transAxes, fontsize=10)
        axFI.text(0.05, 0.15, " ".join(['VgateSlope:', metas[counter]["Status"]["vgateslope"]]), transform=axFI.transAxes, fontsize=10)


        #   sets y label on the left plot
    axFI.set_ylabel('Frequency [Hz]')
    axFreq.set_yticklabels([], visible=False)

    #   sets y lim on both plots
    axFI.set_ylim(axFreq.get_ylim())

    #   writes titles, only over the top figures
    axFI.set_title(subfolder + ": " + "FI rates", fontsize=12)
    axFreq.set_title(subfolder + ": " + "FI freqs", fontsize=12)

    #   writes labels beside the bottom most plot
    axFI.set_xlabel('Current [nA]')
    axFreq.set_xlabel('time [s]')

    #   define file name
    filename = "".join([str(metas[counter]["ReProIndex"]),'_', 'FI_curve'])

    FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
    FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)

    plt.close()

    return filename


def fi_spikes(path_to_folder, subfolder, info):


    """
    fi_spikes parses the information from ficurve-spikes.dat
    :param path_to_folder: location of the experimental folder
    :param subfolder: name of the experimental folder
    :param info: metas of experiment
    :return fnames_rug: raster plot figure names
    :return fnames_FI: file names of the FI plot

    Outputs: graphic files containing FI plots (rates, PSTH, ....) and raster plots

    """
    #   define the input data
    filename = ppjoin(path_to_folder, subfolder, "ficurve-spikes.dat")
    filename2 = ppjoin(path_to_folder, subfolder, "ficurve-data.dat")

    #   get the file, extract the three data subunits
    relacs_file = load(filename)
    relacs_file2 = load(filename2)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        #   load only metas
        metas2,_, _ = relacs_file2.selectall()
    except:
        return None

    ReProIxList = []
    for i in range(len(metas2)):
        ReProIxList.append(metas2[i]["ReProIndex"])

    #   define list for rug filenames
    fnames_rug = []
    fnames_FI = []

    #   for each iteration of the same RePro
    for i in range(len(ReProIxList)):

        #   print processed RePro iteration
        # print("FI Raster & Return ReProIx", ReProIxList[counter[i]])#metas[ReProIx[i]]["ReProIndex"])
        print("FI Raster & Return ReProIx", ReProIxList[i])

        #   select all iterations of the same repro
        metas, _, datas = relacs_file.select({"ReProIndex": ReProIxList[i]})

        #   sort spikes
        spikeDict, stim_amps = fi_spikes_sort(metas, datas)

        #   computes instantaneous spike frequencies based on spikes
        freqDict, timeDict = fi_instant_freq(spikeDict) #, metas[counter[i]:counter[i+1]])
        freq_continuous, freq_continuous2, timeLine, freqSD_continuous = fi_inst_freq_continuous(spikeDict, freqDict, metas)

        #   computes steady state, peak and rest spike frequency, based on above
        fnames2 = fi_curve2(i, path_to_folder, subfolder, timeLine, freq_continuous, freqSD_continuous, metas2, info)

        #   taus, instFreqFit, fitTimeLine = fi_tau_fit(freq_continuous2, metas)
        fnames = fi_rug_plot(i, path_to_folder, subfolder, spikeDict, freq_continuous, freq_continuous2, timeDict,\
                    timeLine, metas, info)

        fnames_rug.append(fnames)
        fnames_FI.append(fnames2)

    return fnames_rug,fnames_FI