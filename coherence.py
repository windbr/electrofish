__author__ = 'jpresern'

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from pyrelacs.DataClasses import load
from os.path import join
from posixpath import join as ppjoin

def process_wn(path_to_folder, fname, data_length):
    """
    Function loads white noise used in stimulation
    :param path_to_folder: folder containting the experimental subfolders
    :param fname: white noise file name joined with the subfolder path
    :param samp_rate: desired sample rate (to match the recording sample rate)
    :param data_length: number of samples in spike train (must match to length white noise)
    :param WN_dur: duration of the white noise
    :return datasWN: vector containing white noise values only.
    """
    # wn_dur /= 1000 #    convert to seconds
    # wn_N_samples = wn_dur*samp_rate

    #   define the input data
    filenameWN = ppjoin(path_to_folder, fname)

    #   get the file, extract the three data subunits
    relacs_file = load(filenameWN)

    #   if relacs file is empty or too short due to aborted RePro presentation
    #   "try:" provides normal termination of the loop instead of error
    try:
        metasWN, _, datasWN = relacs_file.selectall()

    except:
        return None

    #   check if the number of samples matches the number samples in convolved spike train
    datas_whitenoise = datasWN[0][:,1]
    if datas_whitenoise.shape != data_length:
        timepoints_new = np.linspace(datasWN[0][0,0],datasWN[0][-1,0], data_length)
        datas_interpolated = np.interp(timepoints_new, datasWN[0][:,0], datasWN[0][:,1])

    return datas_interpolated


def gauss_kernel(sigma, FS, duration):
    """
    Creates a Gaussian kernel centered in a vector of given duration.

    Poached from the code of dr Jan Grewe

    :param sigma: the standard deviation of the kernel in seconds
    :param sample_rate: the temporal resolution of the kernel in Hz
    :param duration: the desired duration of the kernel in seconds, (in general at least 4 * sigma)

    :return y: the kernel as numpy array
    """
    l = duration * FS
    x = np.arange(-np.floor(l / 2), np.floor(l / 2)) / FS
    y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-(x ** 2 / (2 * sigma ** 2)))

    return y


def train_convolve(train, sigma, FS, wnDur):
    """
    Convolves the spike train with the Gaussian kernel in frequency domain.

    :param train:
    :param sigma: the standard deviation of the kernel in seconds
    :param FS: sampling frequency
    :param wnDur: whit noise duration

    :return conv_Train: convolved spike train
    :return binSpikeTrain: binary spike train

    """
    #   compute the kernel
    kernel = gauss_kernel(sigma, FS, 8 * sigma)

    #   creates binary spike train
    binSpikeTrain, _ = np.histogram(train, (wnDur * FS), (0, wnDur))

    # conv_Train = np.convolve(binSpikeTrain, kernel, mode="same") * FS
    conv_Train = np.convolve(binSpikeTrain, kernel, mode="same")

    return conv_Train, binSpikeTrain

def cohere_transfere_MI (convolved_train, wNoise, nFFT, FS):
    """
    Function computes coherence, mutual information from coherence, transfer function and cross-power-spectra density
    and power spectra density. Function computes the same for whole signal and only for the steady state portion (_short)
    :param metas: meta data
    :param convolved_train: convolved spike train
    :param wNoise: white noise
    :param nFFT: FFT
    :param FS: sample rate
    :return coh: computed coherence
    :return coh_short:
    :return freq:
    :return MI:
    :return MI_short:
    :return P_csd:
    :return P_csd_short:
    :return P_psd:
    :return P_csd_short:
    """

    #   conversions w. clipping first 2 seconds
    wNoiseShort = wNoise[2*FS:]
    conv_train_short = convolved_train[2*FS:]

    #   compute coherence
    coh, freq = mlab.cohere(wNoise-np.mean(wNoise), convolved_train-np.mean(convolved_train), NFFT=nFFT, Fs=FS)
    coh_short, _ = mlab.cohere(wNoiseShort-np.mean(wNoiseShort), conv_train_short-np.mean(conv_train_short), NFFT=nFFT, Fs=FS)

    #   calculate the mutual information
    MI = (-np.log(1-coh))/np.log(2)
    MI_short = (-np.log(1-coh_short))/np.log(2)

    #   compute csd
    P_csd, _ = mlab.csd(wNoise, convolved_train, NFFT=nFFT, Fs=FS)
    P_csd_short, _ = mlab.csd(wNoiseShort, conv_train_short, NFFT=nFFT, Fs=FS)
    if np.any(np.absolute(P_csd_short) < 0):
        print("negative CSD")

    #   compute psd
    P_psd,_ = mlab.psd(wNoise,NFFT=nFFT, Fs=FS)
    P_psd_short,_ = mlab.psd(wNoiseShort,NFFT=nFFT, Fs=FS)

    #   compute the transfer function
    H = np.absolute(P_csd)/P_psd
    H_short = np.absolute(P_csd_short)/P_psd_short

    return freq, coh, coh_short, H, H_short, MI, MI_short, P_csd, P_csd_short, P_psd, P_psd_short

def compute_avgs(coh, coh_short, H, H_short, MI, MI_short):

    #   compute coherence average
    avgCoh = np.average(coh, 0)
    avgCoh = np.reshape(avgCoh, (1,avgCoh.shape[0]))
    avgCoh_short = np.average(coh_short, 0)
    avgCoh_short = np.reshape(avgCoh_short, (1, avgCoh_short.shape[0]))

    #   calculate the mutual information
    avgMI = np.average(MI,0)
    avgMI = np.reshape(avgMI, (1,avgMI.shape[0]))
    avgMI_short = np.average(MI_short,0)
    avgMI_short = np.reshape(avgMI_short, (1,avgMI_short.shape[0]))
    # mut_inf = (-np.log(1-avgCoh))/np.log(2)
    # mut_inf_short = (-np.log(1-avgCoh_short))/np.log(2)

    #   compute the csd, psd average
    avgH = np.average(H,0)
    avgH = np.reshape(avgH, (1, avgH.shape[0]))
    avgH_short = np.average(H_short,0)
    avgH_short = np.reshape(avgH_short, (1, avgH_short.shape[0]))

    #
    # avgP_csd = np.average(np.absolute(P_csd), 0)
    # avgP_psd = np.average(P_psd,0)
    # avgP_csd_short = np.average(np.absolute(P_csd_short), 0)
    # avgP_psd_short = np.average(P_psd_short, 0)
    #
    # #   compute the averaged transfer function
    # avgH = avgP_csd/avgP_psd
    # avgH = np.reshape(avgH, (1, avgH.shape[0]))
    # avgH_short = avgP_csd_short/avgP_psd_short
    # avgH_short = np.reshape(avgH_short, (1, avgH_short.shape[0]))


    return avgCoh, avgCoh_short, avgH, avgH_short, avgMI, avgMI_short

def figure_handles():
    """
    Defines figure and axis handles
    :return: list of handles
    """
    #   defines the figure
    # FHandles = plt.figure(figsize=(14, 10))
    FHandles = plt.figure(figsize=(11.69, 8.27)) # A4


    #   defines coherence axis, number and position
    axCoh = FHandles.add_axes([0.10, 0.4, 0.22, 0.30])
    axCoh_short = FHandles.add_axes([0.10, 0.05, 0.22, 0.30])

    #   defines the transfer function axis
    axH = FHandles.add_axes([0.4, 0.4, 0.22, 0.30])
    axH_short = FHandles.add_axes([0.4, 0.05, 0.22, 0.30])

    #   defines the MI function axis
    axMI = FHandles.add_axes([0.7, 0.4, 0.22, 0.30])
    axMI_short = FHandles.add_axes([0.7, 0.05, 0.22, 0.30])
    #   defines the raster plot axis
    axRaster = FHandles.add_axes([0.10, 0.80, 0.85, 0.15])

    return [FHandles, axCoh, axCoh_short, axMI, axMI_short, axH, axH_short, axRaster]

def plot_the_lot(handles, freq, coh, coh_short, MI, MI_short, H, H_short, metas, cmap, raster = 'empty', \
                 annotation = False, comparison=False):
    """
    :param handles: list of handles for everything
    :parma freq: frequencies
    :param coh:
    :param coh_short:
    :param MI:
    :param MI_short:
    :param H:
    :param H_short:
    :param raster:
    :param color_map:  color map used
    """

    #   define and design grey color map
    # color_spread = np.linspace(0.35,0.8, len(metas))
    # cmap = [ cm.Greys(x) for x in color_spread ]
    lw = 0.25

    #   ensure that there is enough copies of color map values
    if len(cmap) == 1:
        cmap = cmap*coh.shape[0]

    #   use in the single_analysis/parallel_analysis
    if comparison == False and coh.shape[0] == 1:
        lw = 2
        cmap = [cm.Reds(0.6)]

    #   use in the reports
    if comparison == True and coh.shape[0] == 1:
        lw = 2

    #   how many iterations
    for i in range(coh.shape[0]):

        #   plot coherence
        handles[1].plot(freq, coh[i,:], color = cmap[i], linewidth = lw)
        handles[2].plot(freq, coh_short[i,:], color = cmap[i], linewidth = lw)

        #   plot mutual information
        handles[3].plot(freq, MI[i,:], color = cmap[i], linewidth = lw)
        handles[4].plot(freq, MI_short[i,:], color=cmap[i], linewidth= lw)

        #   plot the transfer function
        handles[5].plot(freq, H[i,:], color = cmap[i], linewidth = lw)
        handles[6].plot(freq, H_short[i,:], color = cmap[i], linewidth = lw)

        #   draw raster plot
        if type(raster) == np.ndarray:
            handles[7].plot(raster[i], np.zeros(len(raster[i]))+i, '|', color=cmap[i], ms= 12)
        if type(raster) == list:
            for l in range(len(raster)):
                handles[7].plot(raster[l], np.zeros(len(raster[l]))+l, '|', color=cmap[i], ms= 12)


        #   plot shaping
        wnFreq = float(metas[0]["Settings"]["Waveform"]["freq"].split("H")[0])  # maximum white noise frequency
        wnAmplitude = metas[0]["Settings"]["Waveform"]["amplitude"]  # Standard deviation of the noise
        wnOffset = metas[0]["Settings"]["Stimulus"]["offset"]  # stimulus offset
        wnDC = round(float(metas[0]["Status"]["Current-1"].split('n')[0]),2)  # DC

        #   set x lim
        handles[1].set_xlim(0, wnFreq)
        handles[1].set_ylim(0, 0.6)
        handles[3].set_xlim(0, wnFreq)
        handles[5].set_xlim(0, wnFreq)
        handles[5].set_ylim(0, 600)

        handles[2].set_xlim(0, wnFreq)
        handles[2].set_ylim(0, 0.6)
        handles[6].set_xlim(0, wnFreq)
        handles[6].set_ylim(0, 600)
        handles[4].set_xlim(0, wnFreq)

        #   add labels
        handles[1].set_ylabel("Coherence")
        handles[1].set_xlabel([], visible=False)
        handles[1].set_xticklabels([], visible=False)
        handles[1].set_title("Coherence")
        handles[2].set_xlabel("Frequency [Hz]")
        handles[2].set_ylabel("Coherence")

        handles[5].set_ylabel("Gain [Hz/nA]")
        handles[5].set_xlabel([], visible=False)
        handles[5].set_xticklabels([], visible=False)
        handles[5].set_title("Transfer function")
        handles[6].set_xlabel("Frequency [Hz]")
        handles[6].set_ylabel("Gain [Hz/nA]")

        handles[7].set_xlabel("Time [ms]")
        handles[7].set_ylabel("Iteration")

        handles[3].set_ylabel("MI [bits]")
        handles[3].set_xlabel([], visible=False)
        handles[3].set_xticklabels([], visible=False)
        handles[3].set_title("Mutual information")
        handles[4].set_ylabel("MI [bits]")
        handles[4].set_xlabel("Frequency [Hz]")

    if annotation == True:
        #   comments on the wall
        handles[1].text(0.05, 0.80, " ".join(['RePro Ix:', str(metas[i]["ReProIndex"])]), transform=handles[1].transAxes,
                    fontsize=10)
        handles[1].text(0.05, 0.70, " ".join(['RePro time:', str(metas[i]["ReProTime"])]), transform=handles[1].transAxes,
                    fontsize=10)
        handles[1].text(0.05, 0.60, " ".join(['DC current:', str(wnDC),'nA']), transform=handles[1].transAxes, fontsize=10)
        handles[1].text(0.05, 0.50, " ".join(['Amp. mode:', metas[i]["Status"]["AmplifierMode"]]),
                    transform=handles[1].transAxes, fontsize=10)
        handles[1].text(0.05, 0.40, " ".join(['WN ampp. (SD): ', wnAmplitude]), transform=handles[1].transAxes, fontsize=10)
        handles[1].text(0.05, 0.30, " ".join(['WN offset: ', wnOffset]), transform=handles[1].transAxes, fontsize=10)

        if "SyncPulse" in metas[i]["Status"].keys():

            metas[0]["Status"]["gvgate"].split('.')[0],' nS'
            handles[1].text(0.55, 0.90, " ".join(['DynClamp:', metas[i]["Status"]["SyncPulse"]]), transform=handles[1].transAxes, fontsize=12)
            handles[1].text(0.55, 0.80, " ".join([r'$g_{leak}$',':', metas[i]["Status"]["g"].split('.')[0],'nS']), transform=handles[1].transAxes, fontsize=12)
            handles[1].text(0.55, 0.70, " ".join([r'$E_{leak}$',':', metas[i]["Status"]["E"].split('.')[0],'mV']), transform=handles[1].transAxes, fontsize=12)
            handles[1].text(0.55, 0.60, " ".join([r'$g_{Vgate}$',':', metas[i]["Status"]["gvgate"].split('.')[0],'nS']), transform=handles[1].transAxes, fontsize=12)
            handles[1].text(0.55, 0.50, " ".join([r'$\tau_{Vgate}$',':', metas[i]["Status"]["vgatetau"].split('.')[0],'ms']), transform=handles[1].transAxes, fontsize=12)
            handles[1].text(0.55, 0.40, " ".join([r'$V_{Vgate50}$',':', metas[i]["Status"]["vgatevmid"].split('.')[0],'mV']), transform=handles[1].transAxes, fontsize=12)
            handles[1].text(0.55, 0.30, " ".join([r'$E_{Vgate}$',':', metas[i]["Status"]["Evgate"].split('.')[0],'mV']), transform=handles[1].transAxes, fontsize=12)
            # axCoh.text(0.55, 0.40, " ".join(['VgateSlope:', metas[i]["Status"]["vgateslope"]]), transform=axCoh.transAxes, fontsize=10)

def wht_noise(path_to_folder, subfolder, info):
    """
    Main engine, opens the spike files.

    :param path_to_folder:  folder with recordings
    :param subfolder:       folder containing the experiment/cell
    :param info:            experimental data from info.dat
    :return fnames:         list of figure names to be used in the html build

    """

    #   define sample rate
    FS = 20000

    #   define the Gauss kernel width (= 1SD)
    sigma = 0.001  # seconds, from Berman & Maler 1998

    #   define the input data
    filename = ppjoin(path_to_folder, subfolder, "stimulus-whitenoise-spikes.dat")

    #   get the file, extract the three data subunits
    relacs_file = load(filename)

    #   extract RePro indices
    try:
        ReProIx = relacs_file.fields[('ReProIndex',)]

    except:
        return None

    #   convert set objet into a list
    ReProIx = list(ReProIx)
    ReProIx.sort()

    #   define empty list containing figure names
    fnames = []

    for ix in ReProIx:

        #   if relacs file is empty or too short due to aborted RePro presentation
        #   "try:" provides normal termination of the loop instead of error
        try:
            metas, _, datas = relacs_file.select({"ReProIndex": ix})

        except:
            return None

        print("ReProIx", ix, "Iterations", len(metas))

        #   determine figure handles
        fig = figure_handles()

        #   FFT is defined as something + 1, due to mlab reasoning
        nFFT = 2048
        FFT = (nFFT/2)+1

        #   prepare empty variables
        coh = np.zeros([len(metas), FFT], )
        coh_short = np.zeros([len(metas), FFT], )
        P_csd = np.zeros([len(metas), FFT], dtype=complex)
        P_csd_short = np.zeros([len(metas), FFT], dtype=complex)
        P_psd = np.zeros([len(metas), FFT])
        P_psd_short = np.zeros([len(metas), FFT])
        H = np.zeros([len(metas), FFT],)# dtype=complex)
        H_short = np.zeros([len(metas), FFT],)# dtype=complex)
        MI = np.zeros([len(metas), FFT], )
        MI_short = np.zeros([len(metas), FFT], )
        #   number of stimulus iterations

        for i in range(0, len(metas)):

            color_spread = np.linspace(0.35,0.8, len(metas))
            cmap = [ cm.Greys(x) for x in color_spread ]

            #   extract meta infos
            wnFname = metas[i]["envelope"]
            wnDur = float(metas[i]["Settings"]["Waveform"]["duration"].split("m")[0])  # duration in miliseconds

            #   conversions
            spikes = np.array(datas[i])

            #   conversions
            wnDur /= 1000   #   conversion to miliseconds
            spikes /= 1000

            print(spikes.shape)
            convolved_Train, _ = train_convolve(spikes, sigma, FS, wnDur)
            print(sum(convolved_Train)/FS)
            wNoise = process_wn(path_to_folder, wnFname, len(convolved_Train))

            #   compute coherence, mutual information, transfer and the power spectra and cross-spectra density
            freq, coh[i,:], coh_short[i,:], H[i,:], H_short[i,:], MI[i,:], MI_short[i,:], \
                P_csd[i,:], P_csd_short[i,:], P_psd[i,:], P_psd_short[i,:] \
                = cohere_transfere_MI (convolved_Train, wNoise, nFFT, FS)

        #   plot coherence, mutual information etc....
        plot_the_lot(fig, freq, coh, coh_short, MI, MI_short, H, H_short, metas, cmap, np.array(datas))

        avgCoh, avgCoh_short, avgH, avgH_short, mut_inf, mut_inf_short = compute_avgs(coh, coh_short, H, H_short, MI, MI_short)

        plot_the_lot(fig, freq, avgCoh, avgCoh_short, mut_inf, mut_inf_short, avgH, avgH_short, metas, cmap = [cm.Reds(0.6)],raster='empty', annotation=True)

        fig[2].text(0.05, 0.95, " ".join(['Species:', info["Subject"]["Species"]]), transform=fig[1].transAxes,
                fontsize=10)
        fig[2].text(0.05, 0.90, " ".join(['ELL Segment:', info["Cell"]["Location"]]), transform=fig[1].transAxes,
                fontsize=10)

        #   define file name
        filename = "".join([str(metas[i]["ReProIndex"]), '_', 'whitenoise'])
        fnames.append(filename)

        fig[0].savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'png'])), transparent=True)
        fig[0].savefig(ppjoin(path_to_folder, subfolder, ".".join([filename, 'svg'])), transparent=True)

        plt.close()

    return fnames
