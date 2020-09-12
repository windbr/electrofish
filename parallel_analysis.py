#!/usr/bin/env python3

__author__ = 'jpresern'

"""
Loops through the subfolders containing recordings, performs analysis and outputs the .html file. Uses multiple processors
"""

from joblib import Parallel, delayed
import os, fnmatch
from utils import folder_walk, file_of_interest, seek_unique
from posixpath import join as ppjoin
from collections import OrderedDict, defaultdict
from read_info_dat import read_info

from html_dump import generate_exp_html, generate_index, generate_rug_html

def fn_process(fns):
    """
    Takes list of file names, splits name into path and file name, sorts everything into apropriate semgents
    :param fns: list of lists containing file names with full path (posix)
    :return segment_pages:  dictionary of links to the experiment index pages (MS: MS_index.html; CLS: CLS_index.html,...)
    """

    #   collect all possible locations from the list, extracts every one only once
    segment = seek_unique([fns[i][2].split('_i')[0] for i in range(len(fns))])

    #   define Ordered dict containing links to the experiments in a segment
    segment_pages = OrderedDict()

    #   sort experiment pages by segment
    for seg in segment:
        segment_pages[seg] = [link
                                for link,folder,file in fns
                                if fnmatch.fnmatch(file, str.join("",(seg,'_index.html')))]
    return segment_pages

def exp_process(fns):
    """
    Takes list of file names, links etc. and looks for experiments done in experimental subfolder
    :param fns: list of lists containing file names with full path (posix)
    :return exp_dict: dictionary using links to experiment pages as keys and performed experiments as values
    """
    #   define Ordered dict containing performed experiments
    exp_dict = OrderedDict()

    #   loop over the indexes
    for link, folder, file in fns:
        if folder == './':
            continue
        else:
            filenames = os.listdir(folder)

        ind_exp = []
        if "info.dat" in filenames:
            ind_exp.append(read_info(wd, folder)[0]["Subject"]["Species"])
        else:
            ind_exp.append('unknown')
        if "membraneresistance-trace.dat" in filenames:
            ind_exp.append("Resist")
        if "ficurve-data.dat" in filenames:
            ind_exp.append("FI")
        if "vicurve-data.dat" in filenames:
            ind_exp.append("VI")
        if "stimulus-whitenoise-spikes.dat" in filenames:
            ind_exp.append("Coherence")
        if "saveevents-Spikes-1.dat" in filenames:
            ind_exp.append("Spont act")
        if "transferfunction-data.dat" in filenames:
            ind_exp.append("Transfer")
        # if ind_exp == []:
        #     continue
        # else:
        exp_dict[folder.split('/')[1]] = ind_exp

    return exp_dict



def data_crunch(wd, expfolder):
    os.chdir(ppjoin(wd,expfolder))
    sn = 'single_analysis.py'
    os.system(sn)
    os.chdir('../')

def main():
    """
    Gets the names of subfolders containing experiments, executes single_analysis.py in each of them in parallel,
    collects results and produces index.html in the folder of the execution

    :return:
    """

(wd,expfolder,_) = folder_walk('./')

#   sort the list of experiments
expfolder.sort()

rug_pages = defaultdict(list)
fish_n_chips = OrderedDict()
exp_pages = OrderedDict()

# r = Parallel(n_jobs=1, verbose = 50)(delayed(data_crunch)(wd, expfolder[i]) for i in range(80,84))
# r = Parallel(n_jobs=1, verbose = 50)(delayed(data_crunch)(wd, expfolder[i]) for i in range(80,len(expfolder)))
r = Parallel(n_jobs=7, verbose = 50)(delayed(data_crunch)(wd, expfolder[i]) for i in range(len(expfolder)))
print(r)


# (_,_,filenames) = next(os.walk(os.ppjoin(wd,expfolder)))

htmls = file_of_interest(wd,'*_index.html')
# for i in range(len(htmls)):
#     print(htmls[i])

#   collect dictionary containing links to the experimental pages
exp_pages = fn_process(htmls)
exp_done = exp_process(htmls)
generate_index(exp_pages, wd, exp_done)


if __name__ == '__main__':
    main()
