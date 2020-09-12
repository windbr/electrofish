import numpy as np
import scipy.signal as sig
import scipy.io as sio
import nix
from IPython import embed
import matplotlib.pylab as plt
from sskernel import sskernel

""""

Solutions" for Day 3 "Spectral analysis of spiking responses of
sensory neurons" of the G-Node Winter course in Neural Data Analysis,
2015

by Jan Grewe, Uni Tuebingen
"""

def gauss_kernel(sigma, sample_rate, duration):
    """
    Creates a Gaussian kernel centered in a vector of given duration.

    :param sigma: the standard deviation of the kernel in seconds
    :param sample_rate: the temporal resolution of the kernel in Hz
    :param duration: the desired duration of the kernel in seconds, (in general at least 4 * sigma)

    :return: the kernel as numpy array
    """
    l = duration * sample_rate
    x = np.arange(-np.floor(l / 2), np.floor(l / 2)) / sample_rate
    y = (1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-(x ** 2 / (2 * sigma ** 2)))
    y /= np.sum(y)
    return y


def load_data(filename):
    f = nix.File.open(filename, nix.FileMode.ReadOnly)
    b = f.blocks[0]
    stim_tag = b.tags["stimulus_strong"]
    data = stim_tag.retrieve_data(0)[:]
    sampling_interval = stim_tag.references[0].dimensions[0].sampling_interval
    time = np.asarray(stim_tag.references[0].dimensions[0].axis(data.shape[0]))
    stimulus = stim_tag.retrieve_feature_data(0)[:]
    f.close()
    return data, stimulus, sampling_interval, time


def psd(trace, segment_length, samplerate, apply_window=False, detrend=False):
    nyquist = samplerate/2.
    f_step = samplerate/segment_length
    f = np.arange(0, nyquist + f_step, f_step)
    noOfSamples = len(trace)
    noOfSegments = int(np.floor(noOfSamples/segment_length))
    powers = np.zeros((len(f),noOfSegments))
    window = np.hanning(segment_length)
  
    for n in range(noOfSegments):
        start	= n * segment_length
        end 	= start + segment_length
        segment	= trace[start:end]
        if detrend:
            segment = segment - np.mean(segment)
        if apply_window:
            segment = segment * window
        power       = np.abs(np.fft.fft(segment, segment_length))**2
        powers[:,n] = power[0:np.floor(len(power)/2.) + 1]
    p = np.mean(powers, axis=1)
    return f, p 


def take_a_look_at_the_raw_data(filename):
    data, stimulus, sampling_interval, time = load_data(filename)

    fig = plt.figure()
    # create a rasterplot
    spike_times, trials  = np.nonzero(data)
    spike_times = spike_times * sampling_interval
    ax = fig.add_subplot(311)
    ax.scatter(spike_times, trials, s=1)
    ax.set_xlim([0, 1.])
    ax.set_xlabel("time [s]")
    ax.set_ylabel("trials")

    # plot the PSTH
    kernel = gauss_kernel(0.003, 1. / sampling_interval, 8 * 0.003)
    conv_responses = np.zeros(data.shape)
    for i in range(data.shape[1]):
        conv_responses[:, i] = np.convolve(data[:,i], kernel, mode="same") * 1./sampling_interval
    psth = np.mean(conv_responses, axis=1)
    std_psth = np.std(conv_responses, axis=1)

    ax = fig.add_subplot(312)
    ax.fill_between(time, psth - std_psth, psth + std_psth, color="dodgerblue", alpha=0.25, label="std")
    ax.plot(time, psth, color="dodgerblue", label="response")
    ax.set_xlim([0, 1.])
    ax.set_xlabel("time [s]")
    ax.set_ylabel("firing rate [Hz]")
    ax.legend(fontsize=9)

    # plot the stimulus
    ax = fig.add_subplot(313)
    ax.plot(time, stimulus, label="stimulus", color="silver")
    ax.set_xlim([0, 1.])
    ax.set_xlabel("time [s]")
    ax.set_ylabel("stimulus intensity")
    ax.legend(fontsize=9)

    fig.tight_layout()
    plt.show()

    # plot the cross correlation function
    # do some normalization stuff first, not important for the shape
    r = (psth - np.mean(psth)) / (np.std(psth) * len(psth))
    s = stimulus - np.mean(stimulus) / (np.std(stimulus) * len(stimulus))
                                                        
    c = sig.fftconvolve(s[::-1], r, mode="full") # one could use plt.xcorr for this but is is so much slower
    lags = np.arange(-len(psth), len(psth)-1) * sampling_interval

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(lags, c, color="dodgerblue", lw=2, label="cross correlation")
    ax.set_xlim([-0.1, 0.1])
    ax.set_ylabel("correlation")
    ax.set_xlabel("lag [s]")
    ax.grid(lw=0.75)
    plt.show()


def apply_kde_for_psth(filename):
    data, stimulus, sampling_interval, time = load_data(filename)
    spike_indices, _ = np.nonzero(data)
    spike_times = time[spike_indices]
    _, _, w_opt = sskernel(spike_times)

    kernel = gauss_kernel(w_opt, 1. / sampling_interval, 8 * w_opt)
    conv_responses = np.zeros(data.shape)
    for i in range(data.shape[1]):
        conv_responses[:, i] = np.convolve(data[:,i], kernel, mode="same") * 1./sampling_interval
    psth = np.mean(conv_responses, axis=1)
    std_psth = np.std(conv_responses, axis=1)
    psds = np.zeros((4097,data.shape[1]))
    for i in range(data.shape[1]):
        f, psds[:,i] = psd(conv_responses[:,1], 8192, 1./sampling_interval, detrend=True, apply_window=True)

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(time, psth, color="dodgerblue", label="response, kernel= " + str(w_opt))
    ax.set_xlabel("time [s]")
    ax.set_ylabel("firing rate [Hz]")
    
    ax = fig.add_subplot(212)
    ax.plot(f, np.mean(psds,1), color="dodgerblue")
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("power")
    ax.set_xlim([0, 100])
    plt.show()


def power_spec_surrogate_data():
    samplerate      = 4096.  # Hz
    duration        = 10.  # seconds
    time            = np.arange(1./samplerate, duration, 1./samplerate)   # time vector
    frequencies     = [0.5, 1, 2, 5, 7, 10, 12, 25, 50, 100, 150, 200, 202, 250]
    amplitudes      = np.ones(len(frequencies)) * 0.75
    segment_lengths = [512, 4096, 32768]

    # combine different sinewaves
    y = np.zeros(len(time)) 
    for f, a in zip(frequencies, amplitudes):
        y = y + a * np.sin(2 * np.pi * f * time)
    # add some noise
    y = y + 2*np.random.randn(len(time))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, y, color="dodgerblue")
    ax.set_title('artificial data')
    ax.set_xlabel('time [s]')
    ax.set_ylabel('intensity [arb. units]');
    plt.show()

    # estimate power spectra with different segment lengths
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for s in segment_lengths[::-1]:
        f, p = psd(y, s, samplerate, apply_window=False, detrend=False)
        ax.plot(f, p, label="seg_length="+str(s))
    ax.set_title('power')
    ax.set_xlabel('time [s]')
    ax.set_ylabel('intensity [arb. units]');
    ax.set_xlim([0, 300])
    ax.legend(fontsize=9)
    ax.set_yscale("log")
    plt.show()


def power_spec_real_data(filename, segment_length, apply_window=False):
    """
    Estimate the power spectrum without applying a kernel to smoothen it.
    This exercise should show that the "normal" binary repsesentation of spikes
    may indicate that the neuronal response has power where it should not have.
    """
    data, stimulus, sampling_interval, time = load_data(filename)
    f, stimulus_power = psd(stimulus, segment_length, 1./sampling_interval, apply_window, False)
    # response power - for first trial only
    f, response_power = psd(data[:,1], segment_length, 1./sampling_interval, apply_window, True)
    
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.plot(f, stimulus_power, label="stimulus power")
    ax.set_xlim([1, 1000])
    ax.set_ylim([1, 20000])
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('power')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=9)

    ax = fig.add_subplot(212)
    ax.plot(f, response_power, label="response power")
    ax.set_xlim([1, 1000])
    #ax.set_ylim([1, 1000])
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('power')
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=9)
    plt.show()


def power_spec_real_data_kernels(filename, segment_length):
    """
    Estimate the power spectrum of the data by applying different kernels
    to make the Dirac-like spikes more smooth.
    This exercise should demonstrate how strongly the application of 
    a certain kernel affects the power spectrum and cuts into the range in 
    which the neuron encodes the stimulus.
    """
    kernel_widths = np.asarray([0.25, 1., 2., 5., 10.]) / 1000 # kernel stds in ms
    data, stimulus, sampling_interval, time = load_data(filename)
     
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for k in kernel_widths:
        kernel = gauss_kernel(k, 1. / sampling_interval, 8 * k)
        psds = []
        f = None
        for i in range(data.shape[1]): # each trial
            r = np.convolve(data[:, i], kernel, mode="same")
            f, p = psd(r, segment_length, 1./sampling_interval, True)
            psds.append(p)

        power_spectrum = np.mean(np.asarray(psds), axis=0)
        ax.semilogy(f, power_spectrum, label="kernel width " + str(k))
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("power")
    ax.legend(fontsize=9)
    ax.set_xlim([1, 1000])
    ax.set_ylim([0.0001, 1000])
    plt.show()


def stimulus_response_coherence(filename, segment_length):
    """
    The stimulus response coherence estimated how linearly
    related stimulus and response are. The coherence varies 
    between 0 (no linear correlation) to 1 (perfect linear 
    correlation). Deviations from 1 are attributed to noise 
    or non-linearities in the system.
    To disect these two effects the expected coherence can be
    estimated by calculating the coherence between the individual
    responses and the average (noise-free) resoponse.
    """
    data, stimulus, sampling_interval, time = load_data(filename)
    nyquist = 1./(sampling_interval * 2.)
    f_step = 1./(sampling_interval * segment_length)
    f = np.arange(0, nyquist + f_step, f_step)
    noOfSamples = data.shape[0]
    noOfSegments = int(np.floor(noOfSamples/segment_length))
    kernel = gauss_kernel(0.001, 1./sampling_interval, 0.01)
    window = np.hanning(segment_length)
    coherence_spectra = np.zeros((segment_length, data.shape[1]), dtype=np.complex_)
    exp_coherence_spectra = np.zeros((segment_length, data.shape[1]), dtype=np.complex_)
    # we will need the psth for the expected coherence 
    psth = np.zeros(data.shape[0])
    for i in range(data.shape[1]):
        psth = psth + np.convolve(data[:,i], kernel, mode='same') * (1./sampling_interval)
        psth = psth/data.shape[1]
    # go and calculate the spectra
    for i in range(data.shape[1]):
        trace = data[:,i]/sampling_interval
        trace = np.convolve(trace, kernel, mode="same")
        f_resp = np.zeros((segment_length, noOfSegments), dtype=np.complex_)
        f_psth = np.zeros((segment_length, noOfSegments), dtype=np.complex_)
        f_stim = np.zeros((segment_length, noOfSegments), dtype=np.complex_)
        for n in range(noOfSegments):
            start	= n * segment_length
            end 	= start + segment_length
            resp_segment = trace[start:end]
            resp_segment = resp_segment - np.mean(resp_segment)
            resp_segment = resp_segment * window
            psth_segment = psth[start:end]
            psth_segment = psth_segment - np.mean(psth_segment)
            psth_segment = psth_segment * window
            stim_segment = stimulus[start:end]
            stim_segment = stim_segment - np.mean(stim_segment)
            stim_segment = stim_segment * window
            
            f_resp[:, n] = np.fft.fft(resp_segment, segment_length)
            f_stim[:, n] = np.fft.fft(stim_segment, segment_length)
            f_psth[:, n] = np.fft.fft(psth_segment, segment_length)

        f_resp_conj = np.conjugate(f_resp) # complex conjugate spectrum of response segments
        f_stim_conj = np.conjugate(f_stim) # complex conjugate spectra of stimulus segments
        f_psth_conj = np.conjugate(f_psth) # complex conjugate spectra of psth segments

        sr_cross_spectrum = np.mean(f_stim_conj * f_resp, axis=1) # cross spectrum S*R
        ss_auto_spectrum  = np.mean(f_stim_conj * f_stim, axis=1) # auto spectrum S*S

        rs_cross_spectrum = np.mean(f_resp_conj * f_stim, axis=1) # cross spectrum R*S
        rr_auto_spectrum  = np.mean(f_resp_conj * f_resp, axis=1) # auto spectrum R*R
    
        pr_cross_spectrum = np.mean(f_psth_conj * f_resp, axis=1) # cross spectrum PSTH*R
        pp_auto_spectrum  = np.mean(f_psth_conj * f_psth, axis=1) # auto spectrum PSTH*PSTH
        rp_cross_spectrum = np.mean(f_resp_conj * f_psth, axis=1) # cross spectrum R*PSTH
    
        coherence_spectra[:, i] = (sr_cross_spectrum * rs_cross_spectrum) / (ss_auto_spectrum * rr_auto_spectrum)
        exp_coherence_spectra[:, i] = (pr_cross_spectrum * rp_cross_spectrum) / (pp_auto_spectrum * rr_auto_spectrum)
            
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(f, np.mean(coherence_spectra[:len(f),:], axis=1), color='dodgerblue', label="r-s coherence")
    ax.plot(f, np.mean(exp_coherence_spectra[:len(f),:], axis=1), color='silver', label="r-r coherence")
    ax.set_xlim([0, 300])
    ax.set_ylim([0, 1])
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('coherence')
    ax.legend(fontsize=9)
    plt.show()
    

def reverse_reconstruction_sta(filename):
    data, stimulus, sampling_interval, time = load_data(filename)
    sta_samples = 0.02 / sampling_interval # 20ms before and after the stimulus
    sta = np.zeros(2*sta_samples)
    sta_time_axis = np.arange(-0.02, 0.02, sampling_interval)

    count = np.sum(data)
    spike_indices,_ = np.nonzero(data)
    for index in spike_indices:
        if (index + sta_samples < len(stimulus)) and (index - sta_samples >= 0):
            sta = sta + stimulus[index-sta_samples:index+sta_samples]
        else:
            count -= 1
    sta = sta / count
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sta_time_axis, sta, color='dodgerblue', label="STA")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("stimulus [a.u.]")
    ax.legend(fontsize=9)
    plt.show()
    
    s_estimated = np.zeros(stimulus.shape)
    for i in range(data.shape[1]):
        s_estimated = s_estimated + np.convolve(data[:,i], sta, mode="same")/data.shape[1]
    
    fig =  plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time, s_estimated, color='dodgerblue', label="estimated", zorder=2)
    ax.plot(time, stimulus, color='silver', label="original", zorder=1)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("stimulus [a.u.]")
    ax.set_xlim([0, 0.5])
    ax.legend(fontsize=9)
    plt.show()
   

def reverse_reconstruction_filter(filename, segment_length):
    data, stimulus, sampling_interval, time = load_data(filename)
    nyquist = 1./(sampling_interval * 2.)
    f_step = 1./(sampling_interval * segment_length)
    f = np.arange(0, nyquist + f_step, f_step)
    noOfSamples = data.shape[0]
    noOfSegments = int(np.floor(noOfSamples/segment_length))
    kernel = gauss_kernel(0.001, 1./sampling_interval, 0.01)
    window = np.hanning(segment_length)
    
    psth = np.zeros(data.shape[0])
    for i in range(data.shape[1]):
        psth = psth + np.convolve(data[:,i], kernel, mode='same') * (1./sampling_interval)
        psth = psth/data.shape[1]
  
    filters = np.zeros((segment_length, data.shape[1]), dtype=np.complex_)
    for i in range(data.shape[1]):
        temp = np.convolve(data[:,i], kernel, mode="same")
        f_resp = np.zeros((segment_length, noOfSegments), dtype=np.complex_)
        f_stim = np.zeros((segment_length, noOfSegments), dtype=np.complex_)
        for n in range(noOfSegments):
            start = n * segment_length
            end = start + segment_length
            r_segment = temp[start:end] - np.mean(temp[start:end])
            s_segment = stimulus[start:end] - np.mean(stimulus[start:end])
            f_resp[:,n] = np.fft.fft(r_segment, segment_length)
            f_stim[:,n] = np.fft.fft(s_segment, segment_length)
        f_resp_conj = np.conjugate(f_resp)
        rs_cross_spectrum = np.mean(f_resp_conj * f_stim, axis=1)
        rr_auto_spectrum  = np.mean(f_resp_conj * f_resp, axis=1)
        filters[:,i]  = rs_cross_spectrum / rr_auto_spectrum
    rev_filter = np.mean(filters, axis=1)
    response_spectrum = np.fft.fft(psth, segment_length)
    s_estimated = np.fft.ifft(rev_filter * response_spectrum);
   
    s_estimated = s_estimated* np.std(stimulus)/np.std(s_estimated)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time[:segment_length], s_estimated, color='dodgerblue', label="estimated", zorder=2)
    ax.plot(time[:segment_length], stimulus[:segment_length], color='silver', label="original", zorder=1)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("stimulus [a.u.]")
    ax.set_xlim([0, 0.5])
    ax.set_ylim([-1., 1.])
    ax.legend(fontsize=9)
    plt.show()


if __name__ == "__main__":
    # take_a_look_at_the_raw_data("data_p-unit.h5")
    # power_spec_surrogate_data()
    # power_spec_real_data("data_p-unit.h5", 4096)
    # power_spec_real_data_kernels("data_p-unit.h5", 4096)
    # stimulus_response_coherence("data_p-unit.h5", 8192)
    # reverse_reconstruction_sta("data_p-unit.h5")
    # reverse_reconstruction_filter("data_p-unit.h5", 8192*2)
    apply_kde_for_psth("data_p-unit.h5")
