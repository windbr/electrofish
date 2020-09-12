#!/usr/bin/env python3

__author__ = 'jpresern'

"""
This script processes the Relacs files in the current experimental directory.
It dumps out .html page belonging to the same experiment

"""

# from os.path import join
from posixpath import join as ppjoin
from os import walk, getcwd
from collections import OrderedDict, defaultdict
from read_info_dat import read_info
from FI import fi_curve, fi_spikes, fi_curve2
from VI import vi_curve
from resistance import resistance
from transfer_curve import transferfunc
from coherence import wht_noise
from spont_act import spont_act

from html_dump import generate_exp_html, generate_index, generate_rug_html

#   get the current working dir, split it into the dir and the path to it
wd = getcwd().split(getcwd().split("/")[-1])[0]
expfolder = getcwd().split("/")[-1]

#   define Ordered dictionary in which the .html file names and locations (as keys) are stored
exp_pages = defaultdict(list)
rug_pages = defaultdict(list)


filenames = []
html_page  = []

#   prepare empty dictionary for figure filenames
figs = OrderedDict()
fig_rugs = OrderedDict()

str(expfolder)
#   loads all files
(_,_,filenames) = next(walk(ppjoin(wd, expfolder)))

if "membraneresistance-trace.dat" or "ficurve-data.dat" or "vicurve-data.dat" in filenames:

    if "info.dat" in filenames:
        exp_info = dict(read_info(wd, expfolder)[0])
        print(exp_info)
    else:
        exp_info ={"Cell":{"Location":"UNLABELED"}}

    if "membraneresistance-trace.dat" in filenames:
        print("----------------------->*resistance*")
        figs["resistance"] = resistance(wd, expfolder, exp_info)
        if figs["resistance"] == None:
            del figs["resistance"]
        print("Resistance draw & save done")

    if "ficurve-data.dat" in filenames:
        print("----------------------->*fi*")

        fig_rugs["fi_rug"], figs["fi_curve"] = fi_spikes(wd, expfolder, exp_info)
        if figs["fi_curve"] == None:
             del figs["fi_curve"]

        if fig_rugs["fi_rug"] == None:
            del fig_rugs["fi_rug"]
        else:
            generate_rug_html(wd, expfolder, exp_info, fig_rugs)
        print("FI draw & save done")

    if "vicurve-data.dat" in filenames:
        print("----------------------->*vi*")
        figs["vi_curve"] = vi_curve(wd, expfolder, exp_info)
        if figs["vi_curve"] == None:
            del figs["vi_curve"]
        print("VI draw & save done")

    if "stimulus-whitenoise-spikes.dat" in filenames:
        print("----------------------->*white noise*")
        figs["white_noise"] = wht_noise(wd, expfolder, exp_info)
        if figs["white_noise"] == None:
            del figs["white_noise"]
        print("white noise draw & save done")

    if "saveevents-Spikes-1.dat" in filenames:
        print("---------------------->*spontaneous activity")
        figs["spont_act"] = spont_act(wd, expfolder, exp_info)
        if figs["spont_act"] == None:
            del figs["spont_act"]
        print("spontaneous activity draw & save done")

    if "transferfunction-data.dat" in filenames:
        print("----------------------->*transfer*")
        figs["transfer_function"] = transferfunc(wd, expfolder, exp_info)
        if figs["transfer_function"] == None:
            del figs["transfer_function"]
        print("transfer function draw & save done")


    # collects created page and the segment it belongs to
    html_page, segment = generate_exp_html(wd, expfolder, exp_info, figs)

# else:
#    print('No Relacs data here!')
