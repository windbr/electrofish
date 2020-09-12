#!/usr/bin/env python3

__author__ = 'jpresern'

import sys, getopt, ast, types
import numpy as np
import matplotlib.pyplot as plt
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from os import walk, getcwd
from collections import OrderedDict, defaultdict
from IPython import embed
from utils import command_interpreter


def FI_plot(ReProIx, wd, expfolder, norm=False):

    """"
    Plots the FI curve into the graphic file (*.png, *.svg)
        Requires:
    :param path_to_folder: location of the experimental folder
    :param subfolder: name of the experimental folder
        Outputs:
            graphic files containing FI plots (rates and PSTH) as .png and .svg
    """
    #   define the input data
    filename = ppjoin(wd, expfolder, "ficurve-data.dat")

    #   load data
    relacs_file=load(filename)

  #   four panel figure
    FHandles = plt.figure(figsize=(15, 10))
    axarrPeak   = FHandles.add_axes([0.05, 0.55, 0.40, 0.40])
    axarrSteady = FHandles.add_axes([0.55, 0.55, 0.40, 0.40])
    axarrAvg    = FHandles.add_axes([0.05, 0.05, 0.40, 0.40])
    axarrBase   = FHandles.add_axes([0.55, 0.05, 0.40, 0.40])

    #   define the colormap
    cmapGray = ["Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
           "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue",
           "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue", "Blue", "DarkRed", "DarkOrange","LimeGreen", "Gold", "Plum", "LightSlateGray", "DodgerBlue"]
    counter = 0

    for ix in ReProIx:

        #   if relacs file is empty or too short due to aborted RePro presentation
        #   "try:" provides normal termination of the loop instead of error
        try:
            metas, _, datas = relacs_file.select({"ReProIndex": ix})

        except:
            return None
        #   convert into np array

        #   extract the data for each FI curve
        Istep = datas[0][:, 0]
        FRavg = datas[0][:, 3]
        FRavgSD = datas[0][:, 4]
        FRbase = datas[0][:, 5]
        FRbaseSD = datas[0][:, 6]
        FRpeak = datas[0][:, 9]
        FRpeakSD = datas[0][:, 10]
        FRsteady = datas[0][:, 12]
        FRsteadySD = datas[0][:, 13]

        #   plot data for each FI curve
        axarrBase.plot(Istep, FRbase, color=cmapGray[counter])
        # axarrBase.fill_between(Istep, FRbase+FRbaseSD, FRbase-FRbaseSD, color=cmapGray[counter], alpha=0.2)
        axarrPeak.plot(Istep, FRpeak, color=cmapGray[counter])
        # axarrPeak.fill_between(Istep, FRpeak+FRpeakSD, FRpeak-FRpeakSD, color=cmapGray[counter], alpha=0.2)
        axarrSteady.plot(Istep, FRsteady, color=cmapGray[counter])
        # axarrSteady.fill_between(Istep, FRsteady+FRsteadySD, FRsteady-FRsteadySD, color=cmapGray, alpha=0.2)
        axarrAvg.plot(Istep, FRavg, color=cmapGray[counter])
        # axarrAvg.fill_between(Istep, FRavg+FRavgSD, FRavg-FRavgSD, color=cmapGray[counter], alpha=0.2)

        #   set limits
        axarrPeak.set_ylim(0, 250)
        axarrSteady.set_ylim(0, 200)
        axarrAvg.set_ylim(0, 200)
        axarrBase.set_ylim(0, 200)

        #   set limits
        axarrPeak.set_xlim(-0.8, 0)
        axarrSteady.set_xlim(-0.8, 0)
        axarrAvg.set_xlim(-0.8, 0)
        axarrBase.set_xlim(-0.8, 0)

        #   provide axis labels
        axarrBase.set_xlabel('Current [nA]')
        axarrBase.set_ylabel('Frequency [Hz]')
        axarrPeak.set_xlabel('Current [nA]')
        axarrPeak.set_ylabel('Frequency [Hz]')
        axarrAvg.set_xlabel('Current [nA]')
        axarrAvg.set_ylabel('Frequency [Hz]')
        axarrSteady.set_xlabel('Current [nA]')
        axarrSteady.set_ylabel('Frequency [Hz]')

        #   provide other labels
        axarrPeak.text(0.10, 0.90, " ".join(['Apteronotus leptorhynchus', expfolder]), transform=axarrPeak.transAxes, fontsize = 10)
        axarrPeak.text(0.10, 0.80, " ".join(['DC current:', metas[0]["Status"]["Current-1"]]), transform=axarrPeak.transAxes, fontsize=10)
        axarrPeak.text(0.05, 0.70, " ".join(['Peak current'] ), transform=axarrPeak.transAxes, fontsize=10)
        axarrAvg.text(0.05, 0.70, " ".join(['Average current'] ), transform=axarrAvg.transAxes, fontsize=10)
        axarrSteady.text(0.05, 0.70, " ".join(['Steady current'] ), transform=axarrSteady.transAxes, fontsize=10)
        axarrBase.text(0.05, 0.70, " ".join(['Baseline current'] ), transform=axarrBase.transAxes, fontsize=10)

        #   provide gVgate values
        axarrBase.text(0.05, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["gvgate"]] ), color=cmapGray[counter], transform=axarrBase.transAxes, fontsize=10)
        axarrPeak.text(0.05, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["gvgate"]] ), color=cmapGray[counter], transform=axarrPeak.transAxes, fontsize=10)
        axarrSteady.text(0.05, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["gvgate"]] ), color=cmapGray[counter], transform=axarrSteady.transAxes, fontsize=10)
        axarrAvg.text(0.05, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["gvgate"]] ), color=cmapGray[counter], transform=axarrAvg.transAxes, fontsize=10)

        #   provide gvtau values
        axarrBase.text(0.55, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["vgatetau"]] ), color=cmapGray[counter], transform=axarrBase.transAxes, fontsize=10)
        axarrPeak.text(0.55, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["vgatetau"]] ), color=cmapGray[counter], transform=axarrPeak.transAxes, fontsize=10)
        axarrSteady.text(0.55, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["vgatetau"]] ), color=cmapGray[counter], transform=axarrSteady.transAxes, fontsize=10)
        axarrAvg.text(0.55, 0.05+0.05*counter, " ".join(['g = ', metas[0]["Status"]["vgatetau"]] ), color=cmapGray[counter], transform=axarrAvg.transAxes, fontsize=10)

        #   counter increase
        counter = counter + 1

    #   Save figures
    FHandles.savefig(ppjoin(".".join([expfolder, "dyn_FI_curves", 'png'])), transparent=True)
    FHandles.savefig(ppjoin(".".join([expfolder, "dyn_FI_curves", 'svg'])), transparent=True)

if __name__ == "__main__":
    """
    Deprecated.
    Use report_FI.py instead
    Plots FI curves as computed by relacs.
    Runs from command line in the experiment subfolder.
    Input argument is a list with the RePro indexes which you want to present in the same analysis.
    example: report_FI_curves_dyn.py -l [15,17,21,23]
    """
    #   get the current working dir, split it into the dir and the path to it
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]
    ReProList = command_interpreter(sys.argv[1:])
    print(ReProList["li"])
    FI_plot(ReProList["li"], wd=wd, expfolder=expfolder, norm=True)
