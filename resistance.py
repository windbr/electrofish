__author__ = 'jpresern'


import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from mathEngine import hist_peak_search, baseline_search
from mathEngine import exp_decay

#   TODO: add input resistance measurement!

def resistance(path_to_folder, subfolder, info):

    """
    Resistance plots membrane resistance data. For each RePro run a plot is made
    :param path_to_folder: path to the folder, containing cell/experiment data subfolders
    :param subfolder: subfolder containing cell/experiment data
    :return: Nothing
    :raise: Nothing
    """

    #   define the input data
    filename4 = ppjoin(path_to_folder, subfolder, "membraneresistance-trace.dat")
    filename5 = ppjoin(path_to_folder, subfolder, "membraneresistance-expfit.dat")

    #   get the file, extract the three data subunits
    relacs_file4 = load(filename4)
    relacs_file5 = load(filename5)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        metas4, _, datas4 = relacs_file4.selectall()
        metas5, _, datas5 = relacs_file5.selectall()

    except:
        return None

    #   count the RePro iterations
    resist = np.array(datas4)
    resist_expfit = np.array(datas5)
    ResistIter = len(datas4)

    #   setup plot dimension, 8" width, each RePro repetition gets 3" height
    # TODO: test figure output with a higher dpi settings (220 for retina display, 300 for print)
    # FHandles = plt.figure(figsize=(8, 4*ResistIter))#, dpi=220)

    #   horizontal plot positioning, h size, v size of individual axis pair
    hist_bot, hist_spacer, trace_bot = 0.1, 0.05, 0.22
    hist_width, trace_width = 0.1, 0.75
    # hist_height, trace_height = 0.75/ResistIter, 0.75/ResistIter
    hist_height, trace_height = 0.75, 0.75

     #   empty list for filename output
    fnames = []
    for i in range(0, ResistIter):
        FHandles = plt.figure(figsize=(10, 3))#, dpi=220)
        #   loads meta data for memresistance-trace and mem.resistance-expfit
        comments4 = metas4[i]
        comments5 = metas5[i]

        #   print RePro iteratons
        print("Resistance ReProIx", str(comments4["ReProIndex"]))

        #   defines voltage distribution axis, number and position
        axarrL = FHandles.add_axes([hist_bot, 0.15, hist_width, hist_height])

        #   histogram calculation for the Vrest distribution
        hist, bins = np.histogram(resist[i][resist[i][:,0]<0,1], 50, density=True)
        bincenters = 0.5*(bins[1:]+bins[:-1])
        axarrL.plot(-hist, bincenters, 'dodgerblue', linewidth=1)

        #   find a V_rest
        memPotCalc = hist_peak_search(-hist, bincenters)
        # print(memPotCalc)

        #   defines resistance axis, number and position
        axarrR = FHandles.add_axes([trace_bot, 0.15, trace_width, trace_height])

        #   draws plot
        axarrR.plot(resist[i][:,0], resist[i][:,1], color='dodgerblue')
        axarrR.fill_between(resist[i][:,0], resist[i][:,1]+resist[i][:,2],
                           resist[i][:,1]-resist[i][:,2], color='dodgerblue', alpha=0.2)

        #   set limits on the right plot
        axarrR.set_xlim(min(resist[i][:, 0]), max(resist[i][:, 0]))
        axarrR.set_ylim(min(resist[i][:, 1]) - abs(min(resist[i][:, 1])*0.05),
                        max(resist[i][:, 1]) + abs(max(resist[i][:, 1])*0.05))

        #   set limits on the left plot
        axarrL.set_ylim(axarrR.get_ylim())
        axarrL.set_xlim(-max(hist), min(hist))

        #   define (pre-)stimulus limits
        preStimLim = [min(resist[i][resist[i][:,0]<0,0]), max(resist[i][resist[i][:,0]<0,0])]

        #   draw line marking resting potential obtained by Relacs
        memPotentialRlx = [float(comments4["Vrest"].split('m')[0]), float(comments4["Vrest"].split('m')[0])]
        axarrR.plot(preStimLim, memPotentialRlx, color='r')
        axarrL.axhline(memPotentialRlx[0], color='r', linewidth=1)

        ##   draw line marking resting potential calculated from trace data
        # memPotCalc = [np.round(np.median(resist[i][resist[i][:,0]<0,1])), np.round(np.median(resist[i][resist[i][:,0]<0,1]))]
        if memPotCalc != []:
            axarrR.plot(preStimLim, [memPotCalc[0], memPotCalc[0]], color='m')
            axarrL.axhline(memPotCalc[0], color='gold', linewidth=1)
            axarrR.text(0.25, 0.75, "".join(['V_rest:', str(memPotCalc[0])]), transform=axarrR.transAxes, fontsize=10, color='gold')

        #   if V_rest fails:
        if memPotCalc == []:
            memPotCalc = [0,0]
            memPotCalc[0] = np.median(resist[i][resist[i][:,0]<0,1])
            axarrR.plot(preStimLim, [memPotCalc[0], memPotCalc[0]], color='m')
            axarrL.axhline(memPotCalc[0], color='gold', linewidth=1)
            axarrR.text(0.25, 0.75, "".join(['V_rest:', str(memPotCalc[0])]), transform=axarrR.transAxes, fontsize=10, color='gold')

        #   plot comments and other annotations
        axarrR.text(0.25, 0.95, "".join(['Species:', info["Subject"]["Species"]]), transform=axarrR.transAxes, fontsize = 10)
        axarrR.text(0.25, 0.90, "".join(['ELL Segment:', info["Cell"]["Location"]]), transform=axarrR.transAxes, fontsize = 10)
        axarrR.text(0.25, 0.85, "".join(['RePro Ix:', str(comments4["ReProIndex"])]), transform=axarrR.transAxes, fontsize=10)
        axarrR.text(0.25, 0.8, "".join(['RePro time:', str(comments4["ReProTime"])]), transform=axarrR.transAxes, fontsize=10)
        #axarrR.text(0.25, 0.75, "".join(['V_rest:', comments4["Vrest"]]), transform=axarrR.transAxes, fontsize=10, color='r')
        axarrR.text(0.25, 0.7, "".join(['DC current:', comments4["Status"]["Current-1"]]), transform=axarrR.transAxes, fontsize=10)
        axarrR.text(0.25, 0.65, "".join(['Amp. mode:', comments4["Status"]["AmplifierMode"]]), transform=axarrR.transAxes, fontsize=10)
        axarrR.text(0.25, 0.60, "".join(['Stimulus:', comments4["Settings"]["Stimulus"]["amplitude"]]), transform=axarrR.transAxes, fontsize=10)
        axarrR.text(0.25, 0.20, "".join(['TauM:', comments5["Taum"]]), transform=axarrR.transAxes, fontsize = 10, color='r')
        axarrR.text(0.25, 0.15, "".join(['Resistance:', comments4["Rss"]]), transform=axarrR.transAxes, fontsize = 10, color = 'r')
        if "SyncPulse" in comments4["Status"].keys():
            axarrR.text(0.80, 0.40, "".join(['DynClamp:', comments4["Status"]["SyncPulse"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.35, "".join(['g:', comments4["Status"]["g"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.30, "".join(['E:', comments4["Status"]["E"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.25, "".join(['gVgate:', comments4["Status"]["gvgate"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.20, "".join(['EVgate:', comments4["Status"]["Evgate"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.15, "".join(['VgateTau:', comments4["Status"]["vgatetau"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.10, "".join(['VgateMid:', comments4["Status"]["vgatevmid"]]), transform=axarrR.transAxes, fontsize=10)
            axarrR.text(0.80, 0.05, "".join(['VgateSlope:', comments4["Status"]["vgateslope"]]), transform=axarrR.transAxes, fontsize=10)

        durStim = [float(comments4["Settings"]["Stimulus"]["duration"].split('m')[0]),
                   float(comments4["Settings"]["Stimulus"]["duration"].split('m')[0])]
        axarrR.axvline(0, color='silver')
        axarrR.axvline(durStim[0], color='silver')

        # #   plot Relacs fitted envelope
        # axarrR.plot(resist_expfit[i][:,0], resist_expfit[i][:,1], color='r')

        #   plot internally fit curve
        axarrR.plot(resist_expfit[i][:,0], resist_expfit[i][:,1], color='r')

        #   fits the exponential decay over the voltage response to the current pulse
        # if memPotCalc != []:
        parameters = exp_decay(resist[i][(durStim[0] > resist[i][:,0]) & (0 < resist[i][:,0]),0],
                             resist[i][(durStim[0] > resist[i][:,0]) & (0 < resist[i][:,0]),1], memPotCalc[0])

        fitVoltage = parameters[0]*np.exp(-parameters[1]*resist[i][(durStim[0] > resist[i][:,0]) & (0 < resist[i][:,0]),0])\
                     +parameters[2]

        axarrR.plot(resist[i][(durStim[0] > resist[i][:,0]) & (0 < resist[i][:,0]),0], fitVoltage, color='gold')
        axarrR.text(0.25, 0.50, "".join(['TauM_computed:', str(np.round(1/parameters[1], 1)), 'ms']), transform=axarrR.transAxes, fontsize = 10, color='gold')

        #   find a V_rest as a peakutils.baseline
        # V_base = baseline_search(resist[i][(durStim[0] > resist[i][:,0]) & (0 < resist[i][:,0]),1], 4)
        # axarrR.plot(resist[i][(durStim[0] > resist[i][:,0]) & (0 < resist[i][:,0]),0], V_base, color='blue')

        # #   writes titles
        axarrR.set_title(subfolder + ": " + "Membrane resistance", fontsize=12)

        #   writes x labels beside the bottom most plot
        axarrR.set_xlabel('Time [ms]')
        axarrL.set_xlabel('Distribution')

        #   writes y labels
        axarrL.set_ylabel('Voltage [mV]')
        axarrL.set_xticklabels([], visible=False)
        axarrR.set_yticklabels([], visible=False)

        #   define file name
        filename = "".join([str(comments4["ReProIndex"]),'_', 'Resistance'])
        fnames.append(filename)

        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
        FHandles.savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)
        plt.close()

    return fnames