#!/usr/bin/env python3

__author__ = 'jpresern'

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from posixpath import join as ppjoin
from os import getcwd, walk
from utils import command_interpreter
from FI import dyn_fi_freq
from read_info_dat import read_info

def draw_freaks(ss, peak, rest, metas, repro, info):

    """
    draw_freaks draws steady state and peak frequencies.
    :param ss: ordered dictionary containing a list under ReProIndex. This list contains 1) nd array with desired frequencies paired with currents, gvgate and vgatetau values.
    :param peak: ordered dictionary containing a list under ReProIndex. This list contains 1) nd array with desired frequencies paired with currents, gvgate and vgatetau values.
    :param metas: dict with meta infos
    """
    cmap = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]

    #   define drawing area
    FHandles = plt.figure(figsize=(12, 5.6))
    # FHandles = plt.figure(figsize=(6.43, 3))  #   for publication
    axPeak   = FHandles.add_axes([0.10, 0.10, 0.35, 0.67])
    axSteady = FHandles.add_axes([0.55, 0.10, 0.35, 0.67])
        #   define and design grey color map
    color_spread = np.linspace(0.35,1, len(repro))
    cmapGrey = [ cm.Greys(x) for x in color_spread ]


    #   define counter
    counter = 0

    for k, v in ss.items():

        steady = ss[k][0]
        pk = peak[k][0]
        #   plot without error bars
        axSteady.plot(steady[:,0], steady[:,1], color=cmapGrey[counter], marker='o')
        axPeak.plot(pk[:,0], pk[:,1], color=cmapGrey[counter], marker='o')
        # #   plot with error bars
        # axSteady.errorbar(steady[:,0], steady[:,1], steady[:,2], color=cmap[counter], marker='o')
        # axPeak.errorbar(pk[:,0], pk[:,1], pk[:,2], color=cmap[counter], marker='o')

        if "SyncPulse" in metas[0]["Status"].keys():
            axSteady.text(0.10, 0.90-counter*0.05, "".join([r'$g_{Vgate}$',' = ', ss[k][1].split('.')[0], ' ns']), transform=axSteady.transAxes, fontsize=10, color=cmapGrey[counter])
            axSteady.text(0.750, 0.20-counter*0.05, "".join([r'$\tau_{Vgate}$',' = ', ss[k][2].split('.')[0], ' ms']), transform=axSteady.transAxes, fontsize=10, color=cmapGrey[counter])
            axPeak.text(0.10, 0.90-counter*0.05, "".join([r'$g_{Vgate}$',' = ', peak[k][1].split('.')[0], ' ns']), transform=axPeak.transAxes, fontsize=10, color=cmapGrey[counter])
            axPeak.text(0.750, 0.20-counter*0.05, "".join([r'$\tau_{Vgate}$',' = ', peak[k][2].split('.')[0], ' ms']), transform=axPeak.transAxes, fontsize=10, color=cmapGrey[counter])

        counter = counter +1

    #   plot lims
    # axSteady.set_xlim(-0.8,0.0)
    # axPeak.set_xlim(-0.8,0.0)
    # axSteady.set_ylim(0,250.0)
    # axPeak.set_ylim(0,250.0)

    #   plot labels
    axSteady.set_xlabel("Current [nA]")
    axPeak.set_xlabel("Current [nA]")
    axSteady.set_ylabel("Frequency [Hz]")
    axPeak.set_ylabel("Frequency [Hz]")

    #   titles
    axPeak.set_title("Peak frequency")
    axSteady.set_title("Steady state frequency")
    #   grid
    # axPeak.grid(b=True)
    # axSteady.grid(b=True)

    #     #   extract ReProIx, sort them and merge them
    # keys  =[i for i in ReProList]
    # keys = sorted(list(map(int, keys)))
    # keys = '.'.join(map(str, keys))
    # #   extract current values
    # current = list(ReProList.values())[0]
    # #   join keys, currents and iterations
    name = '_'.join(map(str,repro))

    #   write title on the figure
    FHandles.suptitle("".join([info["Cell"]["Location"],':', expfolder, "FI_dyn",name]), fontsize=12)

    FHandles.savefig(ppjoin(".".join([expfolder, "FI_dyn",name, 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join([expfolder, "FI_dyn",name, 'svg'])), transparent=True)
    #   dump figures into the certain folder, like a data collector folder
    FHandles.savefig(ppjoin('../overviewFI/', ".".join(["_".join([exp_info["Cell"]["Location"],expfolder,"FI_dyn", name]), 'pdf'])), transparent=True)

if __name__ == "__main__":
    """
    Plots and compares FI curves of different RePros for both steady state and peak frequencies
    Computes its own spiking frequencies. Preferred to report_FI_curves_dyn.py
    Runs from command line in the experiment subfolder: report_FI.py
    Input argument is a list with the RePro indexes which you want to present in the same analysis.
    example: report_FI.py -l [15,17,21,23]
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
    print(ReProList["li"])
    ss, peak, metas, rest = dyn_fi_freq(ReProList["li"], wd=wd, expfolder=expfolder)
    draw_freaks (ss, peak, rest, metas, ReProList["li"], exp_info)


