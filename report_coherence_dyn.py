#!/usr/bin/env python3

__author__ = 'jpresern'

import sys
import numpy as np
import matplotlib.cm as cm
from pyrelacs.DataClasses import load
from posixpath import join as ppjoin
from os import getcwd, walk
from utils import command_interpreter
from coherence import process_wn, train_convolve, cohere_transfere_MI, figure_handles, plot_the_lot, compute_avgs
from collections import OrderedDict
from read_info_dat import read_info


def noise_transfer(tests, wd, expfolder):

    """"
    Plots the transfer and coherence curve into the graphic file (*.png, *.svg)
        Requires:
    :param tests:   dictionary containing experimental conditions as keys and ReProIx as values
    :param wd: location of the experimental folder
    :param expfolder: name of the experimental folder

    Outputs:
            graphic files containing coherence, MI and transfer
    """

    #   define sample rate
    FS = 20000

    #   define the Gauss kernel width (= 1SD)
    sigma = 0.001  # seconds, from Berman & Maler 1998

    #   define the input data
    filename = ppjoin(wd, expfolder, "stimulus-whitenoise-spikes.dat")

    #   read in some experimental data
    (_,_,filenames) = next(walk(ppjoin(wd, expfolder)))
    if "info.dat" in filenames:
        exp_info = dict(read_info(wd, expfolder)[0])
        print(exp_info)
    else:
        exp_info ={"Cell":{"Location":"UNLABELED"}}

    #   load data
    relacs_file=load(filename)

    #   four panel figure
    FHandles = figure_handles()

    #   define colors
    color_spread = np.linspace(0.5,0.9, len(tests))
    cmap = [ cm.Greys(x) for x in color_spread ]
    # cmap = [ cm.viridis(x) for x in color_spread ]

    #   color counter
    col_count = 0

     #   FFT is defined as something + 1, due to mlab reasoning
    nFFT = 1024
    FFT = (nFFT/2)+1

    #   define spike dict for the raster plot
    spike_dict  = OrderedDict()

    for k, v in tests.items():
        #   get RePro indexes of the same experimental condition
        ReProIx = tests[k]

        #   define empty variables, containing traces of the same experiments, different RePros
        coh_repro = np.zeros([len(ReProIx), FFT]); coh_repro_short = np.zeros([len(ReProIx), FFT])
        H_repro = np.zeros([len(ReProIx), FFT]); H_repro_short = np.zeros([len(ReProIx), FFT])
        MI_repro = np.zeros([len(ReProIx), FFT]); MI_repro_short = np.zeros([len(ReProIx), FFT])

        spike_list = []

        #   define/reset counter
        counter = 0

        #   iteration of same experimental condition
        for ix in ReProIx:

            try:
                metas, _, datas = relacs_file.select({"ReProIndex": ix})

            except:
                return None

            #   define empty variables
            coh = np.zeros([len(metas), FFT]); coh_short = np.zeros([len(metas), FFT])
            P_csd = np.zeros([len(metas), FFT], dtype=complex); P_csd_short = np.zeros([len(metas), FFT], dtype=complex)
            P_psd = np.zeros([len(metas), FFT]); P_psd_short = np.zeros([len(metas), FFT])
            H = np.zeros([len(metas), FFT]); H_short = np.zeros([len(metas), FFT])
            MI = np.zeros([len(metas), FFT]); MI_short = np.zeros([len(metas), FFT])

            meta_repros = []

            #   number of stimulus iterations
            for i in range(0, len(metas)):

                #   extract meta infos
                wnFname = metas[i]["envelope"]
                wnDur = float(metas[i]["Settings"]["Waveform"]["duration"].split("m")[0])  # duration in miliseconds

                #   spikes
                spikes = np.array(datas[i])

                #   conversions
                wnDur /= 1000   #   conversion to miliseconds
                spikes /= 1000

                print(spikes.shape)
                convolved_Train, _ = train_convolve(spikes, sigma, FS, wnDur)
                print(sum(convolved_Train)/FS)
                wNoise = process_wn(wd, wnFname, len(convolved_Train))

                #   compute coherence, mutual information, transfer and the power spectra and cross-spectra density
                freq, coh[i,:], coh_short[i,:], H[i,:], H_short[i,:], MI[i,:], MI_short[i,:], \
                    P_csd[i,:], P_csd_short[i,:], P_psd[i,:], P_psd_short[i,:] \
                    = cohere_transfere_MI (convolved_Train, wNoise, nFFT, FS)

            #   compute averages over iterations of the *same* repro
            coh_repro[counter,:], coh_repro_short[counter,:], \
                H_repro[counter,:], H_repro_short[counter,:], \
                MI_repro[counter,:], MI_repro_short[counter,:] = compute_avgs(coh, coh_short, H, H_short, MI, MI_short)

            #   store one of the metas
            meta_repros.append(metas[0])
            #   store all the spikes from the same type of experiment
            spike_list.append(datas)
            counter = counter + 1
        #   plot the lot
        plot_the_lot(FHandles, freq, coh_repro, coh_repro_short, MI_repro, MI_repro_short, H_repro, H_repro_short, metas, cmap = [cmap[col_count]], raster='empty', annotation=False, comparison=True)

        #   compute the average of the different repro presentations (with same test conditions)
        avgCoh, avgCoh_short, avgH, avgH_short, avgMI, avgMI_short = compute_avgs(coh_repro, coh_repro_short, H_repro, H_repro_short, MI_repro, MI_repro_short)

        #   plot the lot
        plot_the_lot(FHandles, freq, avgCoh, avgCoh_short, avgMI, avgMI_short, avgH, avgH_short, meta_repros, cmap = [cmap[col_count]], raster='empty', annotation=False, comparison=True)

        #   provide gVgate values
        FHandles[1].text(0.5, 0.8-0.05*col_count, " ".join([r'$g_{Vgate}$',' = ', meta_repros[0]["Status"]["gvgate"].split(".")[0],'nS'] ), color=cmap[col_count], transform=FHandles[1].transAxes, fontsize=10)
        FHandles[5].text(0.5, 0.8-0.05*col_count, " ".join([r'$\tau_{Vgate}$',' = ', meta_repros[0]["Status"]["vgatetau"].split(".")[0],'ms'] ), color=cmap[col_count], transform=FHandles[5].transAxes, fontsize=10)

        #   update spike dictionary
        spike_dict[k] = spike_list

        #   update the color counter
        col_count += 1

    #   write FFT value
    FHandles[1].text(0.05, 0.90, " ".join(['FFT = ',str(nFFT)]), color='k', transform=FHandles[1].transAxes, fontsize=10)

    #   plot raster plot
    spike_iter_count = 0
    for i, k in enumerate(spike_dict):
        for j in range(len(spike_dict[k])):
            for gnj in range(len(spike_dict[k][j])):
                FHandles[7].plot(spike_dict[k][j][gnj], np.zeros(len(spike_dict[k][j][gnj]))+spike_iter_count, '|', color=cmap[i], ms= 12)
                spike_iter_count += 1

    FHandles[7].set_title(ppjoin(".".join([exp_info["Cell"]["Location"],':', expfolder, "dyn_noise_transfer", '_'.join([k for k,v in tests.items()]),
                                              '_'.join([str(x) for sublist in [v for k,v in tests.items()] for x in sublist]),"fft",str(nFFT)])))
    #   Save figures
    FHandles[0].savefig(ppjoin(".".join([expfolder, "coherence_transfer", '.'.join([k for k,v in tests.items()]),
                                        '.'.join([str(x) for sublist in [v for k,v in tests.items()] for x in sublist]),"fft",str(nFFT),'svg'])), transparent=True)
    FHandles[0].savefig(ppjoin(".".join([expfolder, "coherence_transfer", '.'.join([k for k,v in tests.items()]),
                                        ".".join([str(x) for sublist in [v for k,v in tests.items()] for x in sublist]),"fft",str(nFFT),'png'])), transparent=True)
    #   save figures into dedicated folder if necessary
    FHandles[0].savefig(ppjoin('../overviewTransfer/', ".".join(["_".join([exp_info["Cell"]["Location"],expfolder,"coherence_transfer",
                                        "".join([k for k,v in tests.items()]),
                                        "_".join([str(x) for sublist in [v for k,v in tests.items()] for x in sublist]),
                                        "fft",str(nFFT)]), 'pdf'])), transparent=True)

if __name__ == "__main__":

    """
    Compares and plots coherence and transfer curves and MI for different RePro runs. Pool together same experiments in
    from different RePro runs by putting them together in the same arbitrary group.
    Runs from command line in the experiment subfolder: report_coherence_dyn.py
    Input argument is a dictionary with experimental conditions as keys and RePro indexes as a value
    Experimental conditions are arbitrary whatevers, used only to group ReProIndexes, which belong together.
    example: report_coherence_dyn.py -o"{'0nS':[56,58],'400nS':[50,64]}"
    """
    #   get the current working dir, split it into the dir and the path to it
    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]
    ReProList = command_interpreter(sys.argv[1:])
    print(ReProList["od"])
    noise_transfer(ReProList["od"], wd=wd, expfolder=expfolder)
