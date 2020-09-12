#!/usr/bin/env python3

__author__ = 'jpresern'

import sys, getopt, ast, types
import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from os import walk, getcwd
from collections import OrderedDict, defaultdict
from utils import command_interpreter

# return (arg_dict["li"])

def fi_voltage_plot(expfolder):

    #   defines the color map
    cmapBlue = ["Blue", "DarkTurquoise","CadetBlue", "DeepSkyBlue", "CornFlowerBlue", "DodgerBlue", "LightSkyBlue", "LightSteelBlue", "DarkCyan"]

    cmap = ["Blue", "DarkOrange","LimeGreen", "Gold", "Plum", "DarkRed", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
            "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]

    cmapGrey = ["Black", "Grey"]

    #   define the input data file name
    fn1 = "2015-10-14-ab-ficurve4-0.03-3-V.dat"
    fn2 = "2015-10-14-ab-ficurve4-0.03-3-vgate.dat"
    fn3 = "2015-10-14-ab-ficurve31-0.03-3-V.dat"
    fn4 = "2015-10-14-ab-ficurve31-0.03-3-vgate.dat"
    #   get the file and load it
    data1 = np.loadtxt(fn1)
    data2 = np.loadtxt(fn2)
    data3 = np.loadtxt(fn3)
    data4 = np.loadtxt(fn4)

    #   define the absolute size of the raster & return plot
    # FHandles = plt.figure(figsize=(14, 12))
    FHandles = plt.figure(figsize=(6.43, 6))

    #   define individual axes
    # ax4 = FHandles.add_axes([0.1, 0.1, 0.5, 0.1675])
    # ax3 = FHandles.add_axes([0.1, 0.27, 0.5, 0.2175])
    # ax2 = FHandles.add_axes([0.1, 0.60, 0.5, 0.1675])
    # ax1 = FHandles.add_axes([0.1, 0.77, 0.5, 0.2175])
    #
    # ax4 = FHandles.add_axes([0.1, 0.60, 0.7, 0.1675])
    # ax3 = FHandles.add_axes([0.1, 0.77, 0.7, 0.2175])
    # ax2 = FHandles.add_axes([0.1, 0.60, 0.7, 0.1675])
    # ax1 = FHandles.add_axes([0.1, 0.77, 0.7, 0.2175])

    ax4 = FHandles.add_axes([0.15, 0.59, 0.7, 0.12])
    ax3 = FHandles.add_axes([0.15, 0.77, 0.7, 0.2175])
    ax2 = FHandles.add_axes([0.15, 0.59, 0.7, 0.12])
    ax1 = FHandles.add_axes([0.15, 0.77, 0.7, 0.2175])

    #   plot
    ax1.plot(data1[:,0]-500,data1[:,1], color="Black")
    ax2.plot(data2[:,0]-500,-data2[:,1], color="Black")
    ax3.plot(data3[:,0]-500,data3[:,1], color="Grey")
    ax4.plot(data4[:,0]-500,-data4[:,1], color="Grey")

    #   text
    ax1.text(0.10, 0.10, " ".join(['gVgate = 0 nS']), transform=ax1.transAxes, fontsize=10, color="Black")
    ax1.text(0.50, 0.10, " ".join(['gVgate = 60 nS']), transform=ax1.transAxes, fontsize=10, color="Grey")

    #   limits
    ax1.set_xlim(0, 400)
    ax2.set_xlim(0, 400)
    ax3.set_xlim(0, 400)
    ax4.set_xlim(0, 400)

    ax1.set_ylim(-80,-20)
    ax3.set_ylim(-80,-20)
    ax2.set_ylim(-0.05,0.5)
    ax4.set_ylim(-0.05,0.5)

    #   grids
    # ax1.grid(b=True)
    # ax2.grid(b=True)
    # ax3.grid(b=True)
    # ax4.grid(b=True)

    #   axis labels
    ax1.set_xlabel([], visible=False)
    ax3.set_xlabel([], visible=False)
    ax1.set_xticklabels([], visible=False)
    ax3.set_xticklabels([], visible=False)
    ax2.set_xlabel("Time [ms]")
    ax4.set_xlabel("Time [ms]")

    ax2.set_ylabel("Voltage-gated\nurrent [nA]")
    ax1.set_ylabel("Voltage [mV]")
    ax4.set_ylabel("Voltage-gated\ncurrent [nA]")
    ax3.set_ylabel("Voltage [mV]")

    ax4.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4])

    #   Save figures
    FHandles.savefig(ppjoin(".".join([expfolder, "dyn_trace_from_txt", 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join([expfolder, "dyn_trace_from_txt", 'svg'])), transparent=True)
    # FHandles.savefig(ppjoin('../../../../Documents/Znanstveni-prispevki/clanek_DynClamp/', ".".join([expfolder,"dyn_trace_from_txt", 'pdf'])), transparent=True)
    FHandles.savefig(ppjoin('../../../../Documents/Znanstveni-prispevki/clanek_DynClamp/', ".".join(["_".join([expfolder,"dyn_trace_from_txt"]), 'pdf'])), transparent=True)

if __name__ == "__main__":

    """

    """
    #   get the current working dir, split it into the dir and the path to it
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]
    fi_voltage_plot (expfolder)



