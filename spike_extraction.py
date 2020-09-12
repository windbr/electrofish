#!/usr/bin/env python3

__author__ = 'jpresern'

import matplotlib.pyplot as plt
import numpy as np
import sys
from pyrelacs.DataClasses import load
from utils import read_stimuli_position
from os import getcwd
from utils import command_interpreter, detect_peaks
from posixpath import join as ppjoin
from collections import OrderedDict


def data_select(relacs_file, ReProIndex):
    """
    Selectively loads the data belonging to the same RePro
    :param ReProIndex: repro index
    :return: data as a list
    """
    try:
        metas, _, datas = relacs_file.select({"ReProIndex": ReProIndex})
    except:
        return None

    return datas

def prepare_file(expfolder, fn="test.dat", waveform = "Whitenoise"):
    """
    Copies the file "transferfunction-data.dat" into "stimulus-whitenoise-spikes.dat"
    :param fn:  file name
    :return:
    """
    #   create new file
    new_file = open(fn,'w')
    new_file.close()

    within_meta = True
    within_key = False
    within_data = False

    #   prepare counters etc
    counter = 0
    ReProIx = 0

    #   open files
    new_file = open(fn, 'a')

    #   start copying and reforming the file
    with open("transferfunction-data.dat") as open_file:
        for line_no, line in enumerate(open_file):
            # print(line)
            if not line.rstrip().lstrip():  #   is line empty?
                if within_meta == False and within_key == False and counter == 1:   # end of the RePro iteration?
                    new_file.write("\n")
                    within_data = False
                if within_meta == False and within_key == False and counter == 0:
                    new_file.write("\n")
                    counter = counter + 1
                if within_meta == True and within_key == False:
                    #   create new subsection in meta
                    new_file.write("#     Waveform:\n")
                    new_file.write("#         waveform     : "), new_file.write(waveform), new_file.write("\n")
                    new_file.write("#         amplitude    : "), new_file.write(amplitude), new_file.write("\n")
                    new_file.write("#         freq         : "), new_file.write(freq), new_file.write("\n")
                    print(freq)
                    new_file.write("#         duration     : "), new_file.write(duration), new_file.write("\n")
                    new_file.write("\n")
                    within_meta = False
                    continue

            elif line.startswith("#"):
                if within_key == False:
                    #   remember RePro index
                    if line.startswith("# ReProIndex"):
                        ReProIx = line.split(':')[-1]
                        within_meta = True
                        print(ReProIx)
                    #   remember the values for the new subsection of metas
                    if "amplitude" in line:
                        amplitude = line.split(":")[-1].rstrip().lstrip()
                    if "fmax" in line:
                        freq = line.split(":")[-1].rstrip().lstrip()
                        print(freq)
                    if "duration" in line:
                        duration = line.split(":")[-1].rstrip().lstrip()
                    if line.startswith("# ReProTime"):
                        new_file.write("# envelope : ")
                        # new_file.write("\t")
                        new_file.write("\"")
                        new_file.write("".join([expfolder,"/", "Whitenoise","".join([ReProIx.rstrip().lstrip(),".dat"])]))
                        new_file.write("\"")
                        new_file.write("\n")
                    #   are we in the key yet?
                    if line.startswith("#Key"):
                        within_key = True
                        key_start = line_no
                        new_file.write(line)
                        continue
                    new_file.write(line)

                if within_key == True:
                    if key_start == line_no - 1:
                        new_file.write("# t")
                        new_file.write("\n")
                        continue
                    elif key_start == line_no - 2:
                        new_file.write("# ms")
                        new_file.write("\n")
                        continue
                    elif key_start == line_no -3:
                        within_key = False
                        continue
            else:
                if within_meta == False and within_key == False and within_data == False:
                    with open("_".join([ReProIx.rstrip().lstrip(),"spikes.txt"])) as spikes:
                        for spike_line in spikes:
                            new_file.write(spike_line)
                        # new_file.write("\n")
                        # new_file.write("\n")
                        within_data = True
                        continue
                elif within_meta == False and within_key == False and within_data == True:
                    continue
    new_file.close()

if __name__ == "__main__":

    """
    Extracts the spikes from traces. One can use the raw trace or normal (txt) trace. Does not work (yet) with multiple
    iterations of the same stimulus
    Input: list with the repro indexes
    example: spike_extraction -l [2,5,9]
    """


    #   prepare files and metas et al.

    filename = "transferfunction-data.dat"
    filename2 = "transferfunction-traces.dat"
    relacs_file = load(filename)
    relacs_file2 = load(filename2)
    metas, _, _ = relacs_file.selectall()

    wd = getcwd().split(getcwd().split("/")[-1])[0]
    expfolder = getcwd().split("/")[-1]

    # prepare_file(expfolder, fn="stimulus-whitenoise-spikes.dat")

    ReProIxList = command_interpreter(sys.argv[1:])
    if ReProIxList:
        ReProIxList = ReProIxList['li']
    else:
        ReProIxList = [metas[i]["ReProIndex"] for i in range(len(metas))]

    #   voltage threshold
    threshold = 20  #   spike size in milivolts

    trace_dict = OrderedDict()
    spike_dict = OrderedDict()
    # for i in range(len(metas)):
    for i in range(len(ReProIxList)):
        print(ReProIxList[i])
        datas = data_select(relacs_file2, ReProIxList[i])
        # plt.plot(datas[0][:,0], datas[0][:,2])
        pixies, _ = detect_peaks(datas[0][:,0], datas[0][:,2], 20, check_func=None, check_conditions=None )
        # spike_dict[ReProIxList[i]] = pixies
        np.savetxt("_".join([str(ReProIxList[i]),"spikes.txt"]), pixies)
        #   subtract mean from current to get white noise only
        datas[0][:,1] = datas[0][:,1] - np.mean(datas[0][:,1])
        #   trace_dict[ReProIxList[i]] = datas[0][:,0:2]
        np.savetxt("_".join([str(ReProIxList[i]),"whitenoise_traces.txt"]), datas[0][:,0:2])
        #   plt.plot(pixies, np.ones(len(pixies)), '|', color='r', ms= 12)

        #   create new file to contain white noise with header
        wht_noise = open("".join(["Whitenoise",str(ReProIxList[i]),".dat"]),'w')
        wht_noise.close()
        wht_noise = open("".join(["Whitenoise",str(ReProIxList[i]),".dat"]),'a')
        wht_noise.write("# waveform       : Whitenoise"), wht_noise.write("\n")
        wht_noise.write("# frequency      : "), wht_noise.write(metas[i]["Settings"]["Stimulus"]["fmax"]), wht_noise.write("\n")
        wht_noise.write("# amplitude      : "), wht_noise.write(metas[i]["Settings"]["Stimulus"]["amplitude"]), wht_noise.write("\n")
        wht_noise.write("\n")
        wht_noise.write("#Key"), wht_noise.write("\n")
        wht_noise.write("# t   x")
        wht_noise.write("\n")
        wht_noise.write("# s   1")
        wht_noise.write("\n")
        # wht_noise = open("_".join([str(ReProIxList[i]),"whitenoise_traces.dat"]),'a')
        with open("_".join([str(ReProIxList[i]),"whitenoise_traces.txt"])) as open_file:
            for line_no, line in enumerate(open_file):
                wht_noise.write(line)
        wht_noise.close()

    prepare_file(expfolder, fn="stimulus-whitenoise-spikes.dat")

