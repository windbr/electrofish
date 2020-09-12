__author__ = 'jpresern'


import sys, getopt, ast, re
import numpy as np
from scipy.fftpack import hilbert
from collections import OrderedDict
from posixpath import join as ppjoin

def folder_walk(mypath):

    """

    The module manages the folders, containing the experimental data and looks for the relevant files inside these folder.

    folder_walk (mypath)
        returns tuples containing: 1) dir name in question (matching my path), 2) folder names in that
        folder and 3) file names beside them.
        requires: destination path, relative to the folder in which the module resides

    files_of_interest (subfolder)
        returns a list of files inside the subfolder
        requires: subfolder location (join(dirpath,dirnames[x]))
    """

    from os import walk

    (dirpath, dirnames, filenames) = next(walk(mypath))

    return dirpath, dirnames, filenames


def file_of_interest (wd,leftover):

    """
    Browses the subfolders for files matching criteria
    :param wd: directory containing experimental subfolders
    :param sep: separator where to break the file name
    :param leftover: string compared to
    :return: list of file names, matching the criteria
    """

    # subfolder must contain joined path and the subfolder name
    from os import walk
    from posixpath import join as ppjoin
    import fnmatch

    filenames = [[ppjoin(root, name), root, name]
                 for root, dirs, files in walk(wd)
                 for name in files
                 if fnmatch.fnmatch(name, leftover)]


    return filenames


def seek_unique(seq, idfun=None):

    """
    Returns a list of unique values (appearing only once).
    Poached from http://www.peterbe.com/plog/uniqifiers-benchmark
    :param seq: a sequence that needs to be parsed
    :returns result: a list of unique values
    """

    # order preserving
    if idfun is None:
       def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       # if marker in seen: continue
       # seen[marker] = 1
       # result.append(item)
       # improvement by Anonymous in the same web page
       if marker not in seen:
           seen[marker] = 1
           result.append(item)

    return result

def command_interpreter(argv):

    """
    It takes in the commands at the command line and extracts lists, tuples etc.
    py.py --l='[1,2,3,[1,2,3]]' -d "{1:'one',2:'two',3:'three'}" --tu='(1,2,3)'
    :param argv:
    :return:
    """

    arg_dict={}
    switches={'li':list,'di':dict,'tu':tuple,'od':OrderedDict}
    singles=''.join([x[0]+':' for x in switches])
    long_form=[x+'=' for x in switches]
    d={x[0]+':':'--'+x for x in switches}
    try:
        opts, args = getopt.getopt(argv, singles, long_form)
    except getopt.GetoptError:
        print ("bad arg")
        sys.exit(2)

    for opt, arg in opts:
        if opt[1]+':' in d:
            o=d[opt[1]+':'][2:]
        elif opt in d.values():
            o=opt[2:]
        else: o =''

        print (opt, arg, o)
        if o and arg:
            if o == 'od':
                #   dict as it comes
                arg_unordered = ast.literal_eval(arg)

                #   defines ordered keys
                arg_ordered = OrderedDict()

                #   hacks to extract the keys
                arg_key = re.findall(r"'(.*?)':", arg[1:-1], re.DOTALL)

                #   loop over the keys
                for key in arg_key:
                    arg_ordered[key] = arg_unordered[key]

                #   bang key into final structure
                arg_dict[o]=arg_ordered

            else:
                arg_dict[o]=ast.literal_eval(arg)

        if not o or not isinstance(arg_dict[o], switches[o]):
            print (opt, arg, " Error: bad arg")
            sys.exit(2)
    #
    for e in arg_dict:
        print (e, arg_dict[e], type(arg_dict[e]))

    # print(arg_dict["li"])
    # return (arg_dict["li"])
    return (arg_dict)

def envelope(data):
    """
    Envelope of a function. From: Source code for obspy.signal.filter

    Computes the envelope of the given function. The envelope is determined by
    adding the squared amplitudes of the function and it's Hilbert-Transform
    and then taking the square-root. (See [Kanasewich1981]_)
    The envelope at the start/end should not be taken too seriously.

    :type data: numpy.ndarray
    :param data: Data to make envelope of.
    :return: Envelope of input data.
    """
    hilb = hilbert(data)
    data = (data ** 2 + hilb ** 2) ** 0.5
    return data

def detect_peaks( time, data, threshold, check_func=None, check_conditions=None ):
    """
    Function for spike detection
    :param time: data
    :param data: data
    :param threshold: spike size
    :param check_func:
    :param check_conditions:
    :return:
    """


    if not check_conditions:
        check_conditions = dict()

    event_list = list()
    value_list = list()

    # initialize:
    dir = 0
    min_inx = 0
    max_inx = 0
    min_value = data[0]
    max_value = min_value
    trough_inx = 0

    # loop through the new read data
    for index, value in enumerate(data):

        # rising?
        if dir > 0:
            # if the new value is bigger than the old maximum: set it as new maximum
            if max_value < value:
                max_inx = index  # maximum element
                max_value = value

            # otherwise, if the maximum value is bigger than the new value plus the threshold:
            # this is a local maximum!
            elif max_value >= value + threshold:
                # there was a peak:
                event_inx = max_inx

                # check and update event with this magic function
                if check_func:
                    r = check_func( time, data, event_inx, index, trough_inx, min_inx, threshold, check_conditions )
                    if len( r ) > 0 :
                        # this really is an event:
                        event_list.append( r )
                else:
                    # this really is an event:
                    event_list.append( time[event_inx] )
                    #   extract the data value as well
                    value_list.append( data[event_inx] )

                # change direction:
                min_inx = index  # minimum element
                min_value = value
                dir = -1

        # falling?
        elif dir < 0:
            if value < min_value:
                min_inx = index  # minimum element
                min_value = value
                trough_inx = index

            elif value >= min_value + threshold:
                # there was a trough:
                # change direction:
                max_inx = index  # maximum element
                max_value = value
                dir = 1

        # don't know!
        else:
            if max_value >= value + threshold:
                dir = -1  # falling
            elif value >= min_value + threshold:
                dir = 1  # rising

            if max_value < value:
                max_inx = index  # maximum element
                max_value = value

            elif value < min_value:
                min_inx = index  # minimum element
                min_value = value
                trough_inx = index

    return np.array( event_list ), np.array( value_list )

def read_stimuli_position(path_to_folder, subfolder, stim_name, channels = 2):

    """
    Extracts the RePro raw data position and other things from the stimuli.dat
    :param path_to_folder : path to the subfolders, containg experimental data
    :param subfolder : subfolder containing experiment/cell data
    :param stim_name : type of stimulus (RePro) used
    :param ReProIx : list of RePro indices
    :param channels : how many channels to acquire (V, I, VgatedI, etc..) default = 2
    :return raw_data_dict: returns ordered dict where keys are RePro indices and values are lists containing RePro iteration and np.array with values

    """

    filename = ppjoin(path_to_folder, subfolder, "stimuli.dat")
    file = open(filename, "r")
    file_content = file.read()

    #   set counter
    ReProCount = 0

    #   prepare OrderedDict to contain RePro indices as keys
    raw_data_dict = OrderedDict()
    for counter, line in enumerate(file_content.splitlines()):
        if line.find(stim_name) > -1 and line.find("(dataset)") > -1:

            #   data block starts at position + 26
            position = counter + 26
            raw_ix = np.zeros(channels, dtype=float)
            while len(file_content.splitlines()[position]) > 0:
                raw_ix = np.vstack((raw_ix, [np.float(element) for element in file_content.splitlines()[position].split()[0:channels]]))
                position = position + 1

            raw_ix = np.delete(raw_ix, 0, 0)

            raw_data_dict[str(ReProCount)] = raw_ix
            ReProCount = ReProCount + 1

    return(raw_data_dict)





