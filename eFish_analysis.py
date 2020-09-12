#!/usr/bin/env python3

__author__ = 'jpresern'

"""
Loops through the subfolders containing recordings, performs analysis and outputs the .html file. Deprecated.
Use parallel_analysis.py instead
"""
from utils import folder_walk
from os.path import join
from posixpath import join as ppjoin
from os import walk
from collections import OrderedDict, defaultdict
from read_info_dat import read_info
from FI import fi_curve, fi_spikes
from VI import vi_curve
from resistance import resistance
from transfer_curve import transferfunc
from coherence import wht_noise
from spont_act import spont_act

from html_dump import generate_exp_html, generate_index, generate_rug_html
# from ipython import embed

# (wd,expfolder,_) = folder_walk('./')
(wd,expfolder,_) = folder_walk('../Recordings')
# (_,_,filenames) = next(walk(join(wd,expfolder[i])))

#   sort the list of experiments
expfolder.sort()

#   define Ordered dictionary in which the .html file names and locations (as keys) are stored
# exp_pages = defaultdict(list)
rug_pages = defaultdict(list)
fish_n_chips = OrderedDict()
exp_pages = OrderedDict()

# TODO: add fish name to the index.html
# TODO: add RePro (accomplished experiments) names to the index.html

# for i in range(0, len(expfolder)):
for i in range(103, len(expfolder)):
#     i = 43
    filenames = []
    html_page  = []

    #   prepare empty dictionary for figure filenames
    figs = OrderedDict()
    fig_rugs = OrderedDict()

    #   prepare empty dict entry (a list) for species name, experimental RePros used
    fish_n_chips[expfolder[i]] = []

    print(i)
    print(str(expfolder[i]))
    #   loads all files
    (_,_,filenames) = next(walk(ppjoin(wd, expfolder[i])))

    if "membraneresistance-trace.dat" or "ficurve-data.dat" or "vicurve-data.dat" or "transferfunction-data.dat" \
            or "stimulus-whitenoise-spikes.dat" or "saveevents-Spikes-1.dat"in filenames:

        #   extracting the experimental data if available, else replacing them with UNKNOWN
        if "info.dat" in filenames:
            exp_info = dict(read_info(wd, expfolder[i])[0])
            print(exp_info)
            fish_n_chips[expfolder[i]] = [exp_info["Subject"]["Species"]]
        else:
            exp_info ={"Cell":{"Location":"UNLABELED"}, "Subject":{"Species": "UNKNOWN"}}
            print('No Info.dat file, Relacs probably crashed before the experiment was completed!')
            fish_n_chips[expfolder[i]] = [exp_info["Subject"]["Species"]]

        if "membraneresistance-trace.dat" in filenames:
            print("----------------------->*resistance*")
            figs["resistance"] = resistance(wd, expfolder[i], exp_info)
            if figs["resistance"] == None:
                del figs["resistance"]
            print("Resistance draw & save done")
            fish_n_chips[expfolder[i]].extend(["Resistance"])

        if "ficurve-data.dat" in filenames:
            print("----------------------->*fi*")
            fig_rugs["fi_rug"], figs["fi_curve"] = fi_spikes(wd, expfolder[i], exp_info)
            if figs["fi_curve"] == None:
                 del figs["fi_curve"]
            if fig_rugs["fi_rug"] == None:
                del fig_rugs["fi_rug"]
            else:
                generate_rug_html(wd, expfolder[i], exp_info, fig_rugs)
            print("FI draw & save done")
            fish_n_chips[expfolder[i]].extend(["FI"])

        # if "vicurve-data.dat" in filenames:
        #     print("----------------------->*vi*")
        #     figs["vi_curve"] = vi_curve(wd, expfolder[i], exp_info)
        #     if figs["vi_curve"] == None:
        #         del figs["vi_curve"]
        #     print("VI draw & save done")
        #     fish_n_chips[expfolder[i]].extend(["VI"])


        if "stimulus-whitenoise-spikes.dat" in filenames:
            print("----------------------->*white noise*")
            figs["white_noise"] = wht_noise(wd, expfolder[i], exp_info)
            if figs["white_noise"] == None:
                del figs["white_noise"]
            print("white noise draw & save done")
            fish_n_chips[expfolder[i]].extend(["WhiteNoise"])

        if "saveevents-Spikes-1.dat" in filenames:
            print("---------------------->*spontaneous activity")
            figs["spont_act"] = spont_act(wd, expfolder[i], exp_info)
            if figs["spont_act"] == None:
                del figs["spont_act"]
            print("spontaneous activity draw & save done")
            fish_n_chips[expfolder[i]].extend(["spont_act"])

        if "transferfunction-data.dat" in filenames:
            print("----------------------->*transfer*")
            figs["transfer_function"] = transferfunc(wd, expfolder[i], exp_info)
            if figs["transfer_function"] == None:
                del figs["transfer_function"]
            print("transfer function draw & save done")
            fish_n_chips[expfolder[i]].extend(["TransferFunc"])

        # collects created page and the segment it belongs to
        html_page, segment = generate_exp_html(wd, expfolder[i], exp_info, figs)

        # dumps the segment and the page link into the dictionary
        if segment in exp_pages:
            exp_pages[segment].append(ppjoin (".",expfolder[i],html_page))
        else:
            exp_pages[segment] = []
            exp_pages[segment].append(ppjoin (".",expfolder[i],html_page))

generate_index(exp_pages, wd, fish_n_chips)
