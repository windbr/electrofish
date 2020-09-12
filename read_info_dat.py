__author__ = 'jpresern'

from posixpath import join as ppjoin
from pyrelacs.DataClasses import load

def read_info(path_to_folder, subfolder):
    """
    Loads info.dat containing experiment metadata and extracts and returns information of interest
    :param path_to_folder : path to the subfolders, containg experimental data
    :param subfolder : subfolder containing experiment/cell data
    :return meta_info :  returns name of the segment in which the recording was performed

    """

    filename = ppjoin(path_to_folder, subfolder, "info.dat")

    meta_info = load(filename)
    # print(meta_info)
    return(meta_info)