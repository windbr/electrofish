__author__ = 'jpresern'

"""
Digests the transfercurves made by Relacs and plots them.
"""
import numpy as np
import matplotlib.pyplot as plt
# import seaborn as sns
from pyrelacs.DataClasses import load
from os.path import join
from posixpath import join as ppjoin
from itertools import groupby

def transferfunc(path_to_folder, subfolder, info):

    """
    Takes the file called "transferfunction-data.dat" and plots the values found there on the frequency axis. The printed data are "gain (+- SD), phase (+- SD) and coherence (+-SD).
    :param path_to_folder : folder where the experiment subfolders are stored
    :param subfolder : folder containing the experiment itself
    :param info : Data like cell location (map) etc.
    :return fnames : file name containing the figure names
    """

    #   define the input data
    filename1 = ppjoin(path_to_folder, subfolder, "transferfunction-data.dat")
    filename2 = ppjoin(path_to_folder, subfolder, "transferfunction-traces.dat")

    #   get the file, extract the three data subunits
    relacs_file1 = load(filename1)
    relacs_file2 = load(filename2)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        metas1, _, datas1 = relacs_file1.selectall()
        metas2, _, _ = relacs_file2.selectall()
    except:
        return None

    pass


    ReProIxList = []
    for i in range(0,len(metas1)):
        comments1 = metas1[i]
        # extract the unique values only one time
        ReProIxList.append(comments1["ReProIndex"])

    # count occurences of same ReProIx
    ReProIx = [len(list(group)) for key, group in groupby(ReProIxList)]
    ReProIx.insert(0,0) # insert zero for easier job
    # print(ReProIx)

    #   count the iterations
    transfIter = len(datas1)
    print(transfIter)

     #   empty list for filename output
    fnames = []
    counter = 0
    for i in range(0,transfIter):

        # TODO: test figure output with a higher dpi settings (220 for retina display, 300 for print)
        FHandles = plt.figure(figsize=(10, 8))#, dpi=220)

        #   counter increase
        counter = counter + ReProIx[i]
        #   convert the loaded data structures into an appropriate format
        rates1 = np.array(datas1[i])
        datas1[i] = None
        comments1 = metas1[i]
        comments2 = metas2[i]

        #   Frequency limit (only frequency band used in the experiment
        freqLim = float(comments2["Settings"]["Stimulus"]["fmax"].rstrip("Hz"))

        #   TF values - select data
        #   select the data columns
        freq = rates1[:, 0]
        freq = freq[freq<freqLim]
        gain = rates1[0:len(freq), 1]
        gainSD = rates1[0:len(freq), 2]
        phase = rates1[0:len(freq), 3]
        phaseSD = rates1[0:len(freq), 4]
        coh = rates1[0:len(freq), 5]
        cohSD = rates1[0:len(freq), 6]

        #   axes definitions
        axarr = FHandles.add_axes([0.1, 0.5, 0.8, 0.40])
        # axarr = FHandles.add_subplot(1, 1, 1)
        axarr2 = axarr.twinx()
        axarrR = FHandles.add_axes([0.1, 0.1, 0.8, 0.40])

        #   plots transfer curves into the axis
        axarr.plot(freq, gain, color='b', label = "gain")
        axarr.fill_between(freq, gain+gainSD, gain-gainSD, color='b', alpha=0.2)
        axarrR.plot(freq, phase, color='g', label = "phase")
        axarrR.fill_between(freq, phase+phaseSD, phase-phaseSD, color='g', alpha=0.2)
        axarr2.plot(freq, coh, color='r', label = "coherence")
        axarr2.fill_between(freq, coh+cohSD, coh-cohSD, color='r', alpha=0.2)

        #   plot comments
        axarr.text(0.25, 0.95, " ".join(['Species:', info["Subject"]["Species"]]), transform=axarr.transAxes, fontsize = 10)
        axarr.text(0.05, 0.90, " ".join(['ELL Segment:', info["Cell"]["Location"]]), transform=axarr.transAxes, fontsize = 10)
        axarr.text(0.05, 0.85, " ".join(['RePro Ix:', str(comments1["ReProIndex"])]), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.80, " ".join(['RePro time:', str(comments1["ReProTime"])]), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.75, " ".join(['Amp. mode:', comments1["Status"]["AmplifierMode"]]), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.70, " ".join(['Offset:', comments1["Settings"]["Stimulus"]["offset"]]), transform=axarr.transAxes, fontsize=10)
        axarr.text(0.05, 0.65, " ".join(['Amplitude:', comments1["Settings"]["Stimulus"]["amplitude"]]), transform=axarr.transAxes, fontsize=10)
        if "SyncPulse" in comments1["Status"].keys():
            axarr.text(0.05, 0.60, " ".join(['DynClamp:', comments1["Status"]["SyncPulse"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.40, " ".join(['g:', comments1["Status"]["g"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.35, " ".join(['E:', comments1["Status"]["E"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.30, " ".join(['gVgate:', comments1["Status"]["gvgate"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.25, " ".join(['EVgate:', comments1["Status"]["Evgate"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.20, " ".join(['VgateTau:', comments1["Status"]["vgatetau"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.15, " ".join(['VgateMid:', comments1["Status"]["vgatevmid"]]), transform=axarr.transAxes, fontsize=10)
            axarr.text(0.40, 0.10, " ".join(['VgateSlope:', comments1["Status"]["vgateslope"]]), transform=axarr.transAxes, fontsize=10)


        freqLim = float(comments2["Settings"]["Stimulus"]["fmax"].rstrip("Hz"))

        #   sets y label on the left plot
        axarr.set_ylabel('gain (mv/nA)', color = "b")
        axarr2.set_ylabel('coherence', color = "r")
        axarrR.set_ylabel('phase', color = "g")

        #   set x limit
        axarr.set_xlim(0, freqLim)

        #   writes titles, only over the top figures
        axarr.set_title(subfolder + ": " + "transfer function", fontsize=12)

        #   writes labels beside the bottom most plot
        axarrR.set_xlabel('Frequency [Hz]')
        axarr.set_xticklabels([], visible=False)
        axarr2.set_xticklabels([], visible=False)

        #   prepares the legend
        axarr.legend(loc=0)
        axarr2.legend(loc=4)

        #   define file name
        filename = "".join([str(comments1["ReProIndex"]),'_', 'transferfunc'])
        fnames.append(filename)

        #   set color for the y axis
        for tl in axarr2.get_yticklabels():
            tl.set_color("r")
        for tl in axarr.get_yticklabels():
            tl.set_color("b")


        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)

        plt.close()

    return fnames