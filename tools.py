# math
import numpy as np
from scipy.fftpack import fft, rfft, fftshift
from scipy.integrate import trapz
from scipy.signal import medfilt, butter, bessel, firwin, firwin2, lfilter, freqz, decimate, periodogram, welch, iirdesign, get_window, firls, find_peaks
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib
import cmath

# plotting
import matplotlib.pyplot as plt
from matplotlib import axes
ax_obj = axes.Axes
import seaborn as sns

def define_fir_lpf(numtap, cutoff, fs, window=None) :
    nyq = 0.5 * fs
    normalized_cutoff = cutoff / nyq
    if window != None:
        b = firwin(numtap, normalized_cutoff, window=window)
    else:
        b = firwin(numtap, normalized_cutoff)
    a = 1
    return b, a

def calculate_jitter(ssb_pn, fbins, carrier, units='log') :
    '''function to calculate rms jitter (time domain expression of phase noise).
    input: 

    ssb_pn_log- single side band phase noise, relative to the carrier.
    expressed in decibels relative to the carrier in a 1 Hz bandwidth. [dBc/Hz]
    *** dBc = 10*log10(P/Pc). binwidth scaling needs to happen before the logarithm and carrier power normalization.
    fbins- linear frequency bins associated with each phase noise value provided
    carrier- linear frequency value associated with the carrier being referenced for the ssb phase noise values.
    units- choose lienar or logarithmic units of phase noise
    output:
    tj_rms- rms value of the jitter, integrated over bandwidth of the phase noise
    '''

    if units == 'log' : # use the logarithmic values
        ssb_pn_log = np.array(ssb_pn)
        fbins = np.array(fbins)
        ssb_pn_lin = np.power(10, ssb_pn_log/10)
        integrated_pn = 10*np.log10(trapz(y=ssb_pn_lin, x=fbins))
        tj_rms = np.sqrt(2*10**(integrated_pn/10))/(2*np.pi*carrier)

    elif units == 'lin' :    # use the linear values
        ssb_pn_lin= np.array(ssb_pn)
        fbins=np.array(fbins)
        integrated_pn = 10*np.log10(trapz(y=ssb_pn_lin, x=fbins))
        tj_rms = np.sqrt(2*10**(integrated_pn/10))/(2*np.pi*carrier)

    return tj_rms

def plot(d, x, axis = None) :
    ''' simple plot function. supply an axis object to add to an already existing plot.
    *** Recommended to plot less than a million points or matplotlib blows up sometimes. ***
    
    input :
    d : n-length 1-dimensional numpy array
    x :  x-axis
    npt_max : max number of array points to plot to avoid matplotlib crash.
    axis : matplotlib axis object for plotting to.
    '''

    npt = len(x)

    if isinstance(axis, ax_obj) :    # if an axis is supplied, plot to it.
        axis.step(x, d[:npt])

    else :    # no axis, make a quick standalone plot.
        plt.step(x, d[:npt])
        plt.show()

        input('press enter to close plot')
        
        plt.close()

def read_binary(infile, outfile, bits=12, fsr=1.35, send='mid', dtype=np.int16, nchannels=2) :
    '''
    this reads a binary file interpreted as series of 16bit integers, as is the case for our ADC's binary codes.
    two arrays of data are returned, for the adc's channel A and channel B.
    input:
    filename- binary file containing channels A and B. The data is organized:
    (CHA[n=0], CHB[n=0]), (CHA[n=1], CHB[n=1]), ..., (CHA[n=N], CHB[n=N])
    bits- number of bits used by the ADC to construct the data
    fsr- full scale range of the ADC in volts peak to peak
    raw- Set to true to return the raw ADC codes instead of converted voltage values
    '''
    
    maxcode = 2**bits
    midpoint = 2**(bits-1)
        
    if send=='raw' :
        data = np.fromfile(infile, dtype=dtype, count=-1, sep="")
        out=[]
        for i in range(nchannels) :
            out.append(data[i::nchannels])
        out = np.array(out)
    if send=='mid' : 
        dat = np.fromfile(infile, dtype=dtype, count=-1, sep="")
        data = (dat.astype(np.int32)-midpoint)*fsr/maxcode
        out=[]
        for i in range(nchannels) :
            out.append(data[i::nchannels])
        out = np.array(out)
    if send=='float':
        data = np.fromfile(infile, dtype=dtype, count=-1, sep="")
        out=[]
        for i in range(nchannels) :
            out.append(data[i::nchannels])
        out = np.array(out) * fsr/maxcode
        
    np.save(outfile, out)

def open_binary(infile, nchannels=2) :
    data = np.load(infile)
    return data

def make_data(outfile, A, fo, fs, jitter_fs, int_time, n_window) :
    ''' create phase noise corrupted signal with appropriate length 
    input:
    outfile- string for .npy file storing the created signal
    A- signal amplitude
    fo- carrier freq
    fs- sampling freq
    jitter_fs- phase jitter in femtoseconds
    int_time- length of integration window in seconds
    n_window- number of integration windows
    '''
    Wo = 2*np.pi*fo
    l = int(fs*int_time)
    npt = n_window*fs
    carrier= Wo(1.0/fs*np.arange(l) + np.random.normal(loc = 0, scale = t_j*1e-15, size=len(n)))
    # array of samples. 3GSPS over a 1ms window (the helicity rate/integration time).
    n = np.linspace(0, n_window*t_flip, int(n_window*t_flip/T_s))
    #n_avg = int((len(n)-t_off)/n_sample_avg) # number of slices of the data to calculate average amplitude.

    argument = w_c*(n + np.random.normal(loc = 0, scale = t_j, size=len(n)))
    carrier = A_c*np.cos(argument)

def get_freq(data, fs, fo_nominal, int_time, ssb_bw_guess=None, npt2n=False, plot=False, xlim=[10,1.5e9]) :

    '''
    '''

    if npt2n :
        N = int(npt2n)
        binwidth = fs/N
        int_time = N/fs
    else :
        N = int(fs*int_time)    # number of samples in integration window.
        binwidth = (1/int_time) # Hz
    
    bins, v = periodogram(x=data[:N], fs=fs, nfft=None, return_onesided=True, scaling='density')
    tone_i = np.argmax(v)

    # v = rfft(x=data[:N], n=None)
    # v = np.abs(v)
    # center = int(fo_nominal*int_time)
    # left = center - int(ssb_bw_guess*int_time)
    # right = center + int(ssb_bw_guess*int_time)
    
    # tone_i = np.argmax(v[left:right]) + left
    # tone_i2 = np.argmax(v)

    tone_f = tone_i*binwidth
    # bins, power = periodogram(x=data[:N], fs=fs, nfft=None, return_onesided=True, scaling='density')
    # tone_f = np.argmax(power)*binwidth
    f_off = tone_f - fo_nominal
    print(len(data[:N]), 'data points')
    print(len(v), 'fft bins')
    print('frequency resolution = ', binwidth,'Hz\n')
    print('calculated frequency:', tone_f, 'Hz\n')
    # print('calculated frequency', tone_i2*binwidth/2)
    print('frequency offset =', f_off)

    if plot:
        fbins = np.linspace(0, fs/2, int(N/2)+1)
        fig_A, ax_A = plt.subplots(1,1)
        if xlim :
            ax_A.set_xlim(xlim[0], xlim[1])
        else :
            ax_A.set_xlim(1e6, 1e10)
        ax_A.set_xscale('log')
        # ax_A.set_yscale('log')
        ax_A.set_xlabel('Frequency (Hertz)')
        ax_A.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        ax_A.set_title('CH A Noise Spectral Density')
        # ax_A.step(fbins, 10*np.log10(2*np.square(v)/(N*binwidth)))
        ax_A.step(fbins, 10*np.log10(v/np.max(v)))
        fig_A.show()
    
    return tone_f, f_off

def xcor_spectrum(ChA, ChB, fo, fs, nbits=12, int_time=.5e-3, n_window=99, plot=False, dual=False, dualTitle="") :
    '''calculate the fourier spectrum of the adc's digitized signal.
    convert spectrum to units of dBc/Hz (decibels normalized to the carrier power in 1Hz bandwidth).
    calculate the jitter contribution in some specified bandwidth relative to the carrier.
    input:
    
    fo- expected carrier frequency
    fs- ADC sampling clock frequency
    int_time- integration time in seconds.
    int_bw- bandwidth for jitter calculation specified in integer multiples of the binwidth. binwidth = 1/(int_time) Hz
    This bandwidth is single sideband, relative to the carrier (specified in frequency offset, not absolute frequency).
    '''

    # the channel codes are converted to volts, unless 'raw=True'
    N = int(fs*int_time)    # number of samples in integration window.

    binwidth = int(1/int_time) # Hz
    # this indexes the bins to get the desired frequency bins for integrating the phase noise
    # index = np.linspace(int(int_bw[0]), int(int_bw[1]), int(int_bw[1]-int_bw[0])+1, dtype=int)

    Saa = np.zeros(int(N/2)) # power spectrum of channel A
    Sbb = np.zeros(int(N/2)) # power spectrum of channel B
    Sba = np.zeros(int(N/2), dtype=np.complex128) # cross correlation spectrum of channels A and B
    
    #determine the number of samples taken per pulse
    valsPerBin = int(fs*int_time)
    #figure out how many windows to calculate over
    if n_window == "auto":
        n_window = int(np.floor(len(ChA)/valsPerBin))
        print("xcor_spectrum n_window={:}".format(n_window))
    
    for i in range(n_window):
        if np.mod(i,10)==0.:
            print("xcor_spectrum window {} is calculating...".format(i))
        start = int(i*N)
        stop = int((i+1)*N)

        # get positive frequencies of FFT, normalize by N
        a = fft(ChA[start:stop])[:int(N/2)]/N
        b = fft(ChB[start:stop])[:int(N/2)]/N
        
        # sum the uncorrelated variances
        Saa += np.square(np.abs(a))
        Sbb += np.square(np.abs(b))
        Sba += b*np.conj(a)

    # divide by the binwidth and the number of spectrums averaged. multiply by 2 for single sided spectrum.
    # This single-sided power spectral density has units of volts^2/Hz
    Saa = 2*Saa/n_window/binwidth
    Sbb = 2*Sbb/n_window/binwidth

    # each cross correlation spectrum needs complex numbers to be averaged
    # because each calculation uses a complex conjugate.
    # wait to convert to a real PSD until the averaging is complete.
    # this spectrum is due to correlated noise sources.
    Sba = 2*np.abs(Sba)/n_window/binwidth

    fbins = np.linspace(0, fs/2, int(N/2))

    # This spectrum is due to only uncorrelated noise sources.
    Sdiff = (Saa + Sbb - 2*Sba)

    # cutoff0 = 13e6
    # if fo < cutoff0 :
    #     b0, a0 = define_bessel_lpf(cutoff=cutoff0, fs=fs, order=3)
    #     ChA = lfilter(b0, a0, ChA)
    #     ChB = lfilter(b0, a0, ChB)
    if plot:

        fig_A, ax_A = plt.subplots(1,1)

        ax_A.set_xlim(1e6, 1e10)
        ax_A.set_xscale('log')
        #ax_A.set_yscale('log')
        ax_A.set_xlabel('Frequency (Hertz)')
        ax_A.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        ax_A.set_title('CH A Noise Spectral Density')
        ax_A.step(fbins, 10*np.log10(Saa/np.max(Saa)))

        fig_B, ax_B = plt.subplots(1,1)

        ax_B.set_xlim(1e6, 1e10)
        ax_B.set_xscale('log')
        #ax_B.set_yscale('log')
        ax_B.set_xlabel('Frequency (Hertz)')
        ax_B.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        ax_B.set_title('CH B Noise Spectral Density')
        ax_B.step(fbins, 10*np.log10(Sbb/np.max(Sbb)))

        fig_C, ax_C = plt.subplots(1,1)

        ax_C.set_xlim(1e6, 1e10)
        ax_C.set_xscale('log')
        #ax_C.set_yscale('log')
        ax_C.set_xlabel('Frequency (Hertz)')
        ax_C.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        ax_C.set_title('A/B Cross Correlation Noise Spectral Density')
        ax_C.step(fbins, 10*np.log10(np.abs(Sba)/np.max(Sba)))

        fig_D, ax_D = plt.subplots(1,1)

        ax_D.set_xlim(1e6, 1e10)
        ax_D.set_xscale('log')
        #ax_D.set_yscale('log')
        ax_D.set_xlabel('Frequency (Hertz)')
        ax_D.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        ax_D.set_title('Phase Noise Spectral Density')
        ax_D.step(fbins, 10*np.log10(Sdiff/np.max(Sdiff)))

        fig_A.show()
        fig_B.show()
        fig_C.show()
        fig_D.show()
        
    if dual:
        fig_E, ax_E = plt.subplots(1,1)
        ax_E.step(fbins, 10*np.log10(Saa/np.max(Saa)), label='CH A Noise Spectral Density')
        ax_E.step(fbins, 10*np.log10(Sbb/np.max(Sbb)), label='CH B Noise Spectral Density')
        ax_E.set_xscale('log')
        ax_E.set_xlabel('Frequency (Hertz)')
        ax_E.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        ax_E.legend(loc=2)
        ax_E.set_title(dualTitle)
        fig_E.show()

    tone_a = np.argmax(Saa)*binwidth
    tone_b = np.argmax(Sbb)*binwidth

    print('Channel A:', tone_a, 'Hz')
    print('Channel B:', tone_b, 'Hz')

    # tj_A = calculate_jitter(ssb_pn=pn_lin_A[index], fbins=bins_A[index], carrier=fo, units='lin')
    # tj_B = calculate_jitter(ssb_pn=pn_lin_B[index], fbins=bins_B[index], carrier=fo, units='lin')
    return tone_a, tone_b
    input('press enter to finish')

def res_plot(A_avg, B_avg, ax_h=None, ax_s=None) :
    
    '''
    make histogram of the ratio of average amplitudes to show the adc resolution
    make a scatter plot of the average amplitudes of the two adc channels.
    '''

    A = A_avg/B_avg
    n=len(A)

    if ax_h == ax_obj:
        bins, vals, pat = ax_h.hist(x=A, bins=None, range=None)
        ax_h.set_xlabel('(A/B)_avg')        
        ax_h.set_ylabel('number of events')
        ax_h.set_title('(A_avg / B_avg) histogram')
    else :
        fig_h, ax_h = plt.subplots(1,1)
        bins, vals, pat = ax_h.hist(x=A, bins=None, range=None)
        ax_h.set_xlabel('samples')
        ax_h.set_ylabel('volts')
        ax_h.set_title('(A_avg / B_avg) histogram')
        fig_h.show()


    if ax_s == ax_obj:
        ax_s.scatter(x=A_avg, y=B_avg, marker='x')
        ax_s.set_xlabel('Channel A')
        ax_s.set_ylabel('Channel B')
        ax_s.set_title('(average amplitude scatter')    
    else :
        fig_s, ax_s = plt.subplots(1,1)
        ax_s.scatter(x=A_avg, y=B_avg, marker='x')
        ax_s.set_xlabel('Channel A')
        ax_s.set_ylabel('Channel B')
        ax_s.set_title('(average amplitude scatter')
        fig_s.show()
        input('press enter to close resolution plots')

    rms = np.std(A)
    print('resolution = ', rms, 'sigma')


def ddc(ChA, ChB, fo, foff, mixoff, lpf_fc, lpf_ntaps, fs, bits, int_time, ncalc, calc_off, phase_time, nch, 
        plot_en, plot_len, plot_win, plot_Fourier=False, filt='FIR', suppress=True, amp_res=False) :
    
    
    ''' 
    this function calculates the resolution of the ADC for measuring 
    the digitally down converted amplitudes (in 1ms windows)

    input:
    ChA, ChB- numpy arrays of data for the two ADC channels
    fo- linear frequency of the analog signal being sampled
    lpf_fc, lpf_ntaps- cutoff and the number of taps for the lowpass filter.
    fs- sampling frequency of the analog to digital converter
    bits- number of bits precision used. adc has 12bit and 14bit modes.
    int_time- length of each integration window for calculating signal amplitudes in seconds.
    ncalc- number of amplitudes (integration windows) to calculate
    calc_off- number of data points to offset amplitude calculation in each window
    phase_time - length of time, in seconds, to average the phase of the I and Q components
    nch- number of channels to do analysis on, either 1 or 2.
    plot_en- set True to enable plotting
    plot_len- number of data points to plot
    plot_win- index of the integration window to plot.

    output:
    a-
    avg
    avg2
    '''        
    
    fo=fo+foff
    
    Wo = 2*np.pi*(fo-np.abs(mixoff))    
    Ts= 1/fs


    binwidth = int(1/int_time)
    
    #default FIR settings (Hamming window)
    if filt=="FIR":
        cutoff = lpf_fc
        nyq = 0.5 * fs
        b0, a0 = define_fir_lpf(numtap=lpf_ntaps, cutoff=lpf_fc, fs=fs)
        # b2, a2 = define_fir_hpf(numtap=15, cutoff=hpf_fc, fs)
        w, h = freqz(b0, a0, worN=3000000)
    
    #FIR w/ flattop window
    if filt=='flattop':
        cutoff = lpf_fc
        nyq = 0.5 * fs
        b0, a0 = define_fir_lpf(numtap=lpf_ntaps, cutoff=lpf_fc, fs=fs, window=filt)
        # b2, a2 = define_fir_hpf(numtap=15, cutoff=hpf_fc, fs)
        w, h = freqz(b0, a0, worN=3000000)
        
    #FIR w/ Dolph-Chebyshev window
    if filt=='chebwin':
        cutoff = lpf_fc
        nyq = 0.5 * fs
        at = 60
        b0, a0 = define_fir_lpf(numtap=lpf_ntaps, cutoff=lpf_fc, fs=fs, window=(filt, at))
        # b2, a2 = define_fir_hpf(numtap=15, cutoff=hpf_fc, fs)
        w, h = freqz(b0, a0, worN=3000000)
    
    #FIR w/ least squares window
    if filt=='leastsq':
        if np.mod(lpf_ntaps, 2)==0:
            lpf_ntaps = lpf_ntaps + 1
        cutoff = lpf_fc
        nyq = 0.5 * fs
        bands = [cutoff, cutoff*10, cutoff*100, nyq]
        a0 = 1
        b0 = firls(lpf_ntaps, bands, [1,0, 0, 0], nyq=nyq)
        # b2, a2 = define_fir_hpf(numtap=15, cutoff=hpf_fc, fs)
        w, h = freqz(b0, a0, worN=3000000)
    
    #default IIR settings (ellip)
    if filt=="IIR":
        cutoff=lpf_fc*2
        nyq = 0.5 * fs
        normalized_pass = cutoff / nyq
        normalized_stop = 100*cutoff / nyq
        b0, a0 = iirdesign(normalized_pass, normalized_stop, .3, 80)
        w, h = freqz(b0, a0, worN=3000000)
    
    #butterworth IIR
    if filt=="butter":
        cutoff=lpf_fc*2
        nyq = 0.5 * fs
        normalized_pass = cutoff / nyq
        normalized_stop = 100*cutoff / nyq
        b0, a0 = iirdesign(normalized_pass, normalized_stop, .3, 80, ftype=filt)
        w, h = freqz(b0, a0, worN=3000000)
        
    #Chebyshev 1 IIR
    if filt=="cheby1":
        cutoff=lpf_fc*2
        nyq = 0.5 * fs
        normalized_pass = cutoff / nyq
        normalized_stop = 100*cutoff / nyq
        b0, a0 = iirdesign(normalized_pass, normalized_stop, .3, 80, ftype=filt)
        w, h = freqz(b0, a0, worN=3000000)
    
    #Chebyshev 2 IIR
    if filt=="cheby2":
        cutoff=lpf_fc*2
        nyq = 0.5 * fs
        normalized_pass = cutoff / nyq
        normalized_stop = 100*cutoff / nyq
        b0, a0 = iirdesign(normalized_pass, normalized_stop, .3, 80, ftype=filt)
        w, h = freqz(b0, a0, worN=3000000)
        
    #Bessel IIR
    if filt=="bessel":
        cutoff=lpf_fc*2
        nyq = 0.5 * fs
        normalized_pass = cutoff / nyq
        normalized_stop = 100*cutoff / nyq
        b0, a0 = iirdesign(normalized_pass, normalized_stop, .3, 80, ftype=filt)
        w, h = freqz(b0, a0, worN=3000000)
    
    # figf, axf = plt.subplots(1,1)
    plt.plot(0.5*fs*w/np.pi, 10*np.log10(np.abs(h)), 'b')
    plt.plot(np.abs(mixoff), 0.5*np.sqrt(2), 'ko')
    plt.axvline(cutoff, color='k')
    if filt=="IIR":
        plt.axvline(normalized_stop*nyq, color='k')
    #plt.xlim(0, 0.5*fs)
    plt.title("Lowpass Filter Frequency Response")
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Voltage Gain (dB)')
    plt.xscale('log')
    #plt.yscale('log')
    plt.grid()
    plt.show()

    # arrays to store average amplitudes
    avg = np.zeros((nch,ncalc))
    avg2 = np.zeros((nch,ncalc))
    #array for storing raw signal amplitude info
    raw_avgs = np.zeros((nch, ncalc))
    
    # number of samples in integration window
    l = int(int_time*fs)
    npt = ncalc*l # averaging time for the phase. set to integer multiple of the integration window.
    
    Ch = np.array([ChA[:npt], ChB[:npt]])
    
    
# =============================================================================
#     rad = Wo/fs*np.arange(l)
#     I_shift= np.sin(rad) 
#     Q_shift = np.cos(rad)
# =============================================================================

    #try making rad length of full signal
    rad = Wo/fs*np.arange(len(Ch[0]))
    I_shift= np.sin(rad) 
    Q_shift = np.cos(rad)

    I=np.zeros(len(Ch[0]))
    Q=np.zeros(len(Ch[0]))

    for k in range(nch):

        for i in range(ncalc) :

            y = (i*l)
            z = (i+1)*l
            
            # multiply each window (of length l) by the sin and cos modulating terms.
# =============================================================================
#             I[y:z] = Ch[k][y:z]*I_shift
#             Q[y:z] = Ch[k][y:z]*Q_shift
# =============================================================================
            I[y:z] = Ch[k][y:z]*I_shift[y:z]
            Q[y:z] = Ch[k][y:z]*Q_shift[y:z]
            
            #find approximate raw signal amplitude over the bin for diagnostics
            peak = np.max(Ch[k][y:z])
            trough = np.min(Ch[k][y:z])
            raw_avgs[k][i] = (peak - trough)/2

        # low pass filter frequency down-mixed data
        I_f = lfilter(b0, a0, I)
        Q_f = lfilter(b0, a0, Q)

        phase_npt = int(phase_time*fs)
        
        for i in range(ncalc) :
        
            d = int(i*l + int(calc_off))  # this offset is to avoid the filter's settling time, which 
            e = int((i+1)*l)              # is on the order of the number of taps for the FIR filter.
            f = int(i*l + phase_npt) 
            
            avg_phis = np.arctan(I_f[d:f]/Q_f[d:f])
            #avg_phi = np.mean(avg_phis)
            
            # phase based reconstruction
            a = np.abs(2*(Q_f[d:e]*np.cos(avg_phis) + I_f[d:e]*np.sin(avg_phis)))
        
            # pythagorean reconstruction
            a2 = 2*np.hypot(I_f[d:e], Q_f[d:e])

            # average amplitude recovered from I,Q components
            avg[k][i] = np.mean(a)
            avg2[k][i] = np.mean(a2)
    
    if plot_en :
        
        if (filt=="IIR" or filt=="butter" or filt=="cheby1" or filt=="cheby2" or filt=="ellip" or filt=="bessel"):
            begin=0
        else:
            begin=lpf_ntaps
        off = plot_win*l
        plot_len = int(plot_len)
        xaxis = np.arange(begin, plot_len)
        plt.rcParams['agg.path.chunksize'] = 10000
        if suppress==False:
            print('sample spacing = {0:.5f}... nanoseconds'.format(Ts*1e9))
        end = off + plot_len
        figA, axA = plt.subplots(1,1)
        axA.set_xlabel('samples')
        axA.set_ylabel('volts')
        axA.set_title('adc channel A raw data')
        axA.plot(xaxis, ChA[begin:end])
        figA.show()

        figB, axB = plt.subplots(1,1)
        axB.set_xlabel('samples')
        axB.set_ylabel('volts')
        axB.set_title('adc channel B raw data')
        axB.plot(xaxis, ChB[begin:end])
        figB.show()

        fig1, ax1 = plt.subplots(1,1)    
        ax1.set_xlabel('samples')
        ax1.set_ylabel('amplitude')
        ax1.set_title('I component')
        ax1.plot(xaxis, I[begin:end])
        fig1.show()

        fig2, ax2 = plt.subplots(1,1)
        ax2.set_xlabel('samples')
        ax2.set_ylabel('amplitude')
        ax2.set_title('Q component')
        ax2.plot(xaxis, Q[begin:end])
        fig2.show()

        fig3, ax3 = plt.subplots(1,1)
        ax3.set_xlabel('samples')
        ax3.set_ylabel('amplitude')
        ax3.set_title('I component (filtered)')
        ax3.plot(xaxis, I_f[begin:end])
        fig3.show()

        fig4, ax4 = plt.subplots(1,1)
        ax4.set_xlabel('samples')
        ax4.set_ylabel('amplitude')
        ax4.set_title('Q component (filtered)')
        ax4.plot(xaxis, Q_f[begin:end])
        fig4.show()

        
        fig5, ax5 = plt.subplots(1,1)
        ax5.set_xlabel('samples')
        ax5.set_ylabel('amplitude')
        ax5.set_title('Phase Reconstructed Amplitude')
        ax5.plot(xaxis[begin:len(a)], a[begin:end])
        fig5.show()
        
        fig6, ax6 = plt.subplots(1,1)
        ax6.set_xlabel('samples')
        ax6.set_ylabel('amplitude')
        ax6.set_title('Hypot Reconstructed Amplitude')
        ax6.plot(xaxis[begin:len(a2)], a2[begin:end])
        fig6.show()
        
    if plot_Fourier:
        
        # the channel codes are converted to volts, unless 'raw=True'
        N = int(fs*int_time)    # number of samples in integration window.
        fbins = np.linspace(0, fs/2, int(N/2))
    
        binwidth = int(1/int_time) # Hz
        # this indexes the bins to the desired frequency bins for integrating the phase noise
        # index = np.linspace(int(int_bw[0]), int(int_bw[1]), int(int_bw[1]-int_bw[0])+1, dtype=int)
        
        Sa = np.zeros(int(N/2)) # power spectrum of A
        Sb = np.zeros(int(N/2)) # power spectrum of B
        SsI = np.zeros(int(N/2)) # power spectrum of sin mixer term
        SsQ = np.zeros(int(N/2)) # power spectrum of cos mixer term
        Si = np.zeros(int(N/2)) # power spectrum of I
        Sq = np.zeros(int(N/2)) # power spectrum of Q
        Sif = np.zeros(int(N/2)) # power spectrum of I_f
        Sqf = np.zeros(int(N/2)) # power spectrum of Q_f
        
        n_window = int(np.floor(len(I)/N))
        
        if suppress==False:
            print('Doing Fourier transforms...')
            
        for i in range(n_window):
            
            if suppress==False:
                print('...')
            
            start = int(i*N)
            stop = int((i+1)*N)
    
            # get positive frequencies of FFT, normalize by N
            af = fft(ChA[start:stop])[:int(N/2)]/N
            bf = fft(ChB[start:stop])[:int(N/2)]/N
            Ssif = fft(I_shift[start:stop])[:int(N/2)]/N
            Ssqf = fft(Q_shift[start:stop])[:int(N/2)]/N
            i = fft(I[start:stop])[:int(N/2)]/N
            q = fft(Q[start:stop])[:int(N/2)]/N
            i_f = fft(I_f[start:stop])[:int(N/2)]/N
            q_f = fft(Q_f[start:stop])[:int(N/2)]/N
            
            # sum the uncorrelated variances
            Sa += np.square(np.abs(af))
            Sb += np.square(np.abs(bf))
            SsI += np.square(np.abs(Ssif))
            SsQ += np.square(np.abs(Ssqf))
            Si += np.square(np.abs(i))
            Sq += np.square(np.abs(q))
            Sif += np.square(np.abs(i_f))
            Sqf += np.square(np.abs(q_f))
    
        # divide by the binwidth and the number of spectrums averaged. multiply by 2 for single sided spectrum.
        # This single-sided power spectral density has units of volts^2/Hz
        Sa = 2*Sa/n_window/binwidth
        Sb = 2*Sb/n_window/binwidth
        SsI = 2*SsI/n_window/binwidth
        SsQ = 2*SsQ/n_window/binwidth
        Si = 2*Si/n_window/binwidth
        Sq = 2*Sq/n_window/binwidth
        Sif = 2*Sif/n_window/binwidth
        Sqf = 2*Sqf/n_window/binwidth
        
        figg, axg = plt.subplots(1,1)
        axg.step(fbins, 10*np.log10(Sa/np.max(Sa)), label='ChA Spectral Density')
        axg.step(fbins, 10*np.log10(Sb/np.max(Sb)), label='ChB Spectral Density')
        axg.step(fbins, 10*np.log10(SsI/np.max(SsI)), label='Mixer sin term Spectral Density')
        axg.step(fbins, 10*np.log10(SsQ/np.max(SsQ)), label='Mixer cos term Spectral Density')
        axg.set_xscale('log')
        axg.set_xlabel('Frequency (Hertz)')
        axg.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        axg.legend(loc=2)
        axg.set_title("Power Spectrums")
        figg.show()
        
        figf, axf = plt.subplots(1,1)
        axf.step(fbins, 10*np.log10(Si/np.max(Si)), label='I Spectral Density')
        axf.step(fbins, 10*np.log10(Sq/np.max(Sq)), label='Q Spectral Density')
        axf.step(fbins, 10*np.log10(Sif/np.max(Sif)), label='I filtered Spectral Density')
        axf.step(fbins, 10*np.log10(Sqf/np.max(Sqf)), label='Q filtered Spectral Density')
        axf.set_xscale('log')
        axf.set_xlabel('Frequency (Hertz)')
        axf.set_ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        axf.legend(loc=3)
        axf.set_title("Power Spectrums")
        figf.show()
    
    #starting at second pulse window calculate amplitude asymmetry resolution of phase reconstructed amplitude
    if amp_res:
        
        A = avg[0][2:]
        B = avg[1][2:]
        
        calcs = len(A)
        
        Apairs = []
        Bpairs = []
        
        for i in range(0, calcs, 2) :
            
            #start with A
            a1 = A[i]
            a2 = A[i+1]
            Apairs.append([a1,a2])
            
            #same process for B
            b1 = B[i]
            b2 = B[i+1]
            Bpairs.append([b1,b2])
        
        #convert lists to arrays to make calculations easier
        Apairs = np.array(Apairs)
        Bpairs = np.array(Bpairs)
    
        #calculate relative differences within channels
        Adiffs = (Apairs[:,0] - Apairs[:,1])/(Apairs[:,0] + Apairs[:,1])
        Bdiffs = (Bpairs[:,0] - Bpairs[:,1])/(Bpairs[:,0] + Bpairs[:,1])
        
        #calculate differences between channels
        diffs = Adiffs - Bdiffs
    
        #stat params
        sigma = np.std(diffs)
        mu = np.mean(diffs)
        
        print('Amplitude asymmetry resolution: {:.3f} ppm'.format(1e6*sigma))
        
        #histogram plot the difference in amplitude asymmetry differences between channels
        figh, axh = plt.subplots(1,1)
        y,x,_ = axh.hist(diffs*1e6, range = (-int(3*sigma*1e6), int(3*sigma*1e6)), bins=20)
        axh.set_xlabel('Amplitude Asymmetry Difference (PPM)')
        axh.text(-int(2.5*sigma*1e6), y.max(), s='Sigma={0:.3f} PPM'.format(sigma*1e6))
        axh.set_title('Phase reconstructed amplitude asymmetry resolution')
        figh.show()
        
        figi, axi = plt.subplots(1,1)
        axi.scatter(Adiffs*1e6, Bdiffs*1e6)
        axi.set_xlabel('Channel A asymmetry (PPM)')
        axi.set_ylabel('Channel B asymmetry (PPM)')
        axi.set_title('Phase reconstructed amplitude asymmetry resolution')
        figi.show()

    return a, avg, avg2, raw_avgs

def calibrate(A, B, fsr=1.35, bits=12, ax1=None, ax2=None) :
    ''' make a 2D scatter plot of A and B.
    Fit a line to it to get the slope.
    input:
    A- channel A
    B- channel B

    output:
    slope- fitted slope of the line
    '''
    # sort the numpy arrays so we can fit them.
    o = np.argsort(A)
    A = A[o]
    B = B[o]

    # 1D polynomial fit (a line)
    z = np.polyfit(A, B, 1)
    fit = np.poly1d(z)
    slope = z[0]
    offset= z[1]
    f = lambda x : offset + slope*x


    print('line has form:', offset, '+', slope, 'x')

    f1, a1 = plt.subplots(1,1)
    a1.scatter(A,B, color='blue', label='A/B')
    a1.plot(A,  f(A), color='red', label='fit')
    a1.set_xlabel('A (volts)')
    a1.set_ylabel('B (volts)')
    a1.set_title('A vs B scatter')
    a1.legend()
    f1.show()

    f2, a2 = plt.subplots(1,1)
    diff = B-(slope*A + offset)
    val, bins, pat = a2.hist(diff, bins=100)
    a2.set_xlabel('A (volts)')
    a2.set_ylabel('A - k*B + Vo (volts)')
    a2.set_title('calibrated A/B difference histogram')
    f2.show()

    # sigma = np.std(val*bins)
    print(bits, 'bit mode.' , fsr, 'volts full scale range.\n')
    print('1 bit precision = ', fsr/2**bits, 'volts')
    print('2 bit precision = ', fsr/2**(bits-1), 'volts')
    print('3 bit precision = ', fsr/2**(bits-2), 'volts\n')
    # print('A-B distribution:', sigma, 'volts sigma')

    input('press enter to finish calibration')
    return val, bins, slope, offset
    
def amplitude_asym_hist(ChA, ChB, fs=3e9, ncalc="auto", pulse_width=5e-4, hist=True, scatter=False, title="", suppress=False, ppm_plot=False, diffs_plot=False):
    
    
    ''' 
    this function calculates the amplitude asymmetry between successive pairs of integration windows.
    the difference in this quantity between channels A and B is calculated over a number of bins and histogram plotted with a gaussian fit

    input:
    ChA, ChB- numpy arrays of data for the two ADC channels
    fs- sampling frequency of the analog to digital converter
    bits- number of bits precision used. adc has 12bit and 14bit modes.
    ncalc- number of pulse widths to use in calculation
    pulse_width - length of time for a single bin. in an experiment would correspond to +/- helicity pulse time. give in seconds
    hist - determines whether or not to plot a histogram of amplitude asymmetries
    scatter - determines whether or not to plot a scatter plot of amplitude asymmetries in channels A & B

    output:
    sigma - standard dev in resolution
    mu - average amplitude asymmetry difference
    diffs - array of successive amplitude asymmetry differences
    '''  

    #determine the number of samples taken per pulse
    valsPerBin = int(fs*pulse_width)

    #lists to store Channel A & B amplitude pairs
    Apairs = []
    Bpairs = []
    
    if ncalc == "auto":
        ncalc = int(np.floor(len(ChA)/valsPerBin))
        if suppress==False:
            print("ncalc={:}".format(ncalc))

    #populate Apairs & Bpairs with bin averaged RMS amplitudes
    for i in range(0, ncalc, 2) :

        #start with A
        a1 = ChA[i*valsPerBin:(i+1)*valsPerBin-1]
        a2 = ChA[(i+1)*valsPerBin:(i+2)*valsPerBin-1]
        a1RMS = np.sqrt(np.dot(a1,a1))
        a2RMS = np.sqrt(np.dot(a2,a2))
        Apairs.append([a1RMS,a2RMS])

        #same process for B
        b1 = ChB[i*valsPerBin:(i+1)*valsPerBin-1]
        b2 = ChB[(i+1)*valsPerBin:(i+2)*valsPerBin-1]
        b1RMS = np.sqrt(np.dot(b1,b1))
        b2RMS = np.sqrt(np.dot(b2,b2))
        Bpairs.append([b1RMS,b2RMS])
        
        if (np.mod(i, 10)==0 and suppress==False):
            print('Calculating window {} in amplitude resolution...'.format(i))

    #convert lists to arrays to make calculations easier
    Apairs = np.array(Apairs)
    Bpairs = np.array(Bpairs)

    #calculate relative differences within channels
    Adiffs = (Apairs[:,0] - Apairs[:,1])/(Apairs[:,0] + Apairs[:,1])
    Bdiffs = (Bpairs[:,0] - Bpairs[:,1])/(Bpairs[:,0] + Bpairs[:,1])
    
    if ppm_plot:
        fppm, appm = plt.subplots(1,1)
        appm.plot(Adiffs*1e6, label="channel A asymmetries")
        appm.plot(Bdiffs*1e6, label="channel B asymmetries")
        appm.set_title(title)
        fppm.show()
    
    #calculate differences between channels
    diffs = Adiffs - Bdiffs
    
    if diffs_plot:
        fppm, appm = plt.subplots(1,1)
        appm.plot(diffs*1e6, label="diffs")
        appm.set_title(title)
        fppm.show()
    
    #stat params
    sigma = np.std(diffs)
    mu = np.mean(diffs)

    #histogram plot the difference in amplitude asymmetry differences between channels
    if hist :
        f1, a1 = plt.subplots(1,1)
        y,x,_ = a1.hist(diffs*1e6, range = (-int(3*sigma*1e6), int(3*sigma*1e6)), bins=20)
        a1.set_xlabel('Amplitude Asymmetry Difference (PPM)')
        a1.text(-int(2.5*sigma*1e6), y.max(), s='Sigma={0:.3f} PPM'.format(sigma*1e6))
        a1.set_title(title)
        f1.show()
    
    #scatter plot
    if scatter :
        f2, a2 = plt.subplots(1,1)
        a2.scatter(Adiffs*1e6, Bdiffs*1e6)
        a2.set_xlabel('Channel A asymmetry (PPM)')
        a2.set_ylabel('Channel B asymmetry (PPM)')
        a2.set_title(title)
        f1.show()

    return sigma, mu, diffs

def plot_res_vs_binwidth(ChA, ChB, minWidth, maxWidth, steps, title, ncalc = 30, fs=3e9, log=False) :
    
    ''' 
    calculates resolution as a function of bin width and plots it

    input:
    ChA, ChB- numpy arrays of data for the two ADC channels
    fs- sampling frequency of the analog to digital converter
    minWidth - minimum pulse_width to calculate resolution for
    maxWidth - maximum pulse_width to calculate resolution for
    steps - number of pulse_widths to calculate resolutions for, using minWidth and maxWidth as range limits
    ncalc - number of pulse pairs to resolution over
    bits- number of bits precision used. adc has 12bit and 14bit modes.
    '''
    
    #instantiate binwidth values and create an array to store resolution data
    widthList = np.logspace(np.log10(minWidth), np.log10(maxWidth), num=steps)
    resList = np.zeros_like(widthList)
    
    #loop to calculate resolutions at various bin widths
    for i in range(0,len(widthList)) :
        if np.mod(i, 300)==0:
            print('Calculating resolution for {0:f} millisecond pulse windows'.format(1e3*widthList[i]))
        res, mu, diffs = amplitude_asym_hist(ChA, ChB, fs=fs, ncalc=ncalc, pulse_width=widthList[i], hist=False, scatter=False, suppress=True)
        resList[i] = res
    
    fig, ax = plt.subplots(1,1)
    if log==False :
        ax.plot(widthList*1e3, resList*1e6, '.')
        ax.set_xlabel('Pulse width (ms)')
        ax.set_ylabel('Resolution (PPM)')
        ax.set_title(title)
        fig.show()
    else:
        ax.plot(widthList*1e3, resList*1e6, '.')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('Pulse width (ms)')
        ax.set_ylabel('Resolution (PPM)')
        ax.set_title(title)
        fig.show()
    
    return resList, widthList

def power_scatter(fileList, adcZ, fs=3e9, bits=12, pulse_width=5e-4, title=""):
    '''
    Takes a list of filevoltage.npy files and the ADC impedance. Returns the calculated power,
    the calculated amplitude asymmetry difference, and creates a scatter plot.
    
    inputs:
    fileList: list of filevoltage.npy filenames to calculate v_rms from
    adcZ: impedance (in Ohms) of ADC input
    fs: sampling frequency
    bits: bit precision
    pulse_width: time per bin
    '''
    
    #will store calculated power & resolution information for each .npy file
    powerResList = []
    
    #iterates through list and calculates the signal power & resolution
    for i in range(0,len(fileList)):
        
        #unpack the file
        A, B = open_binary(fileList[i])
        
        #RMS squared voltage of each channel
        v_rms2A = np.dot(A,A) / len(A)
        v_rms2B = np.dot(B,B) / len(B)
        
        #power for each channel
        powerA = v_rms2A/adcZ
        powerB = v_rms2B/adcZ
        
        powerAve = (powerA+powerB)/2
        
        #convert to dBm
        powerdB = 10 * np.log10(powerAve/1e-3)
        
        #calculate amplitude asymmetry diff
        sigma, mu, diffs = amplitude_asym_hist(A, B, fs=fs, bits=bits, ncalc="auto", pulse_width=pulse_width, hist=False)
        
        #store the newly calculated info in powerResList
        powerResList.append([powerdB, sigma*1e6])
    
    #plot up the data
    plt.scatter(np.array(powerResList)[:,0], np.array(powerResList)[:,1])
    plt.title(title)
    plt.xlabel('Signal Power (dBm)')
    plt.ylabel('Amplitude Asymmetry Resolution (PPM)')
    plt.ylim(bottom=0)
    #plt.yscale('log')
    plt.show()
    
    return powerResList

def DDC_amplitudes(A_I, A_Q, B_I, B_Q, plot=True, title="") :
    """
    Takes four output streams from a DDC signal and recombines them into two streams
    giving the signal amplitudes for channels A&B
    
    Inputs:
    A_I/A_Q/B_I/B_Q: down mixed streams outputted by DDC
    plot: boolean instruction whether to plot amplitudes for first 800 samples
    title: title for plot
    
    Outputs:
    A_amp: hypot reconstructed amplitude of channel A
    B_amp: hypot reconstructed amplitude of channel B
    """
    
    #calculate hypotenuse reconstructed amplitudes pointwise along channels A and B
    A_amp = 2*np.hypot(A_I, A_Q)
    B_amp = 2*np.hypot(B_I, B_Q)
    
    #plot it
    if plot:
        fig, ax = plt.subplots(1,1)
        ax.plot(A_amp[0:800])
        ax.plot(B_amp[0:800])
        ax.set_title(title)
        ax.set_ylabel("Hypotenuse reconstructed amplitude")
        ax.set_xlabel("Sample number")
        fig.show()
    
    return A_amp, B_amp
    
def downmix(signal, fs, ms, single=None, offset=0) :
    """
    Down mixes a signal and returns the I and Q data streams.
    
    Inputs:
    signal: input signal data
    fs: signal sampling frequency
    ms: mixer frequency
    fixAmp: if true, returned values of I and Q are doubled to account for halving of amplitude during downmixing
    single: instead of downmixing to Q and I, can just dowmix to Q or I. If single=Q mix to Q only, if single=I mix to I only
        if single= "complex" returns single, complex downmixed stream
    offset: phase offset
    
    Outputs:
    I: signal * sin
    Q: signal * cos
    """
    
    #length of signal in time
    deltat = len(signal) / fs
    
    #array of time values for producing mixer signal
    time = np.linspace(0, deltat, len(signal))
    
    if single=="I":
        
        #calculate & return I
        mixI = np.cos(2*np.pi*time*ms+offset)
        return signal*mixI
    
    elif single=="Q":
        
        #calculate & return Q
        mixQ = np.sin(2*np.pi*time*ms+offset)
        return signal*mixQ
    
    elif single == "complex" :
        
        #calculate complex mixer signal, mix & return mixed signal
        mix = np.exp(1j*2*np.pi*ms*time)
        return signal * mix
        
    else:
        
        #mixer signals
        mixQ = np.sin(2*np.pi*time*ms+offset)
        mixI = np.cos(2*np.pi*time*ms+offset)
        
        #calculate I & Q
        I = signal * mixI
        Q = signal * mixQ
        
        return I, Q

def single_sideband(signal, fs, ms, dual=False) :
    """
    Takes a downmixed signal and returns the lower and upper sidebands
    
    Inputs:
    signal: input signal data.
        If inputing a downmixed dual stream format as signal=[I, Q] and set dual=True
    fs: downmixed signal sampling frequency
    ms: mixer frequency
    decim: downconversion factor
    
    Outputs;
    Hssb: higher sideband signal
    Lssb: lower sideband signal
    """
    
    if dual:
        I = signal[0]
        Q = signal[1]
        
    else:
        I, Q = downmix(signal, fs, ms)
        
    II, IQ = downmix(I, fs, ms)
    QI, QQ = downmix(Q, fs, ms)
    
    #produce original downmixed signal and its Hilbert transform
    s = II-QQ
    shilb = -(IQ+QI)
    
    #modulation
    Scos, Ssin = downmix(s, fs, ms)
    shilbCos, shilbSin = downmix(shilb, fs, ms)
    
    #calculate upper and lower sidebands
    Hssb = Scos - shilbSin
    Lssb = Scos + shilbSin
    
    return Hssb, Lssb

def FFT_plot(signalList, fs, int_time, scale='log', labels=False, returnFigAx=False) :
    """
    Plots the complex FFT of a array of signals
    
    Inputs:
    signalList: list of signals to plot. Formatted like [signal1, ..., signalN]
    fs: sampling frequency of signals
    int_time: binwidth in time
    labels: determines whether plot will have labels & legend. If not false, formatted like ["label 1", ..., "label n"].
        - Must be same length as signalList
    returnFigAx: if true, doesn't plot, but returns fig and ax for further modification
    """
    
    N = int(fs*int_time)    # number of samples in integration window.

    binwidth = int(1/int_time) # Hz
    
    #determine the number of samples taken per pulse
    valsPerBin = int(fs*int_time)
    
    S = []
    fig, ax = plt.subplots(1,1)
    fbins = np.linspace(-fs/2, fs/2, N)
    for j in range(len(signalList)):
        
        #temporary variable to hold spectrum for jth signal
        Sj =  np.zeros(N)
        
        #figure out how many windows to calculate over
        n_window = int(np.floor(len(signalList[j])/valsPerBin))
        
        print('Working on plotting signal {}'.format(j+1))
        
        for i in range(n_window):
            start = int(i*N)
            stop = int((i+1)*N)
    
            # get frequencies of FFT, normalize by N
            a = fft(signalList[j][start:stop])/N
            
            #shift frequencies into proper bins for plotting
            a = fftshift(a)
            
            # sum the uncorrelated variances
            Sj += np.abs(a)
            
        # divide by the binwidth and the number of spectrums averaged
        # This single-sided power spectral density has units of volts^2/Hz
        Sj = Sj/n_window/binwidth
        Sj = Sj**2 # convert to power spectrum
        if labels!= False:
            ax.step(fbins, 10*np.log10(Sj/np.max(Sj)), label=labels[j])
            
        else:
            ax.step(fbins, 10*np.log10(Sj/np.max(Sj)))
        
    ax.set_xscale(scale)
    ax.set_xlabel('Frequency (Hertz)')
    ax.set_ylabel(r'$\frac{dBm}{Hz}$', rotation=0, fontsize=16)
    if labels!=False:
        ax.legend(loc=2)
    
    if returnFigAx==False:    
        fig.show()
    else:
        return fig, ax
    
def decimate(signal, decim):
    """
    Takes a signal and downsamples it
    
    Inputs:
    signal: input to be downsampled
    decim: decimation factor
    
    Output:
    downsamp: downsampled signal
    """
    
    downsamp = signal[::int(decim)]
    
    return downsamp

# =============================================================================
# def phase_reconstruction(I, Q, fs, fo, ms):
#     """
#     Performs quadrature phase reconstruction on a downmixed & filtered signal
#     
#     Inputs:
#     Q: cos modulated signal
#     I: sin modulated signal
#     fs: sample frequency
#     fo: signal frequency
#     ms: mixer frequency
#     
#     Outputs:
#     amps: array of amplitude at each point
#     """
#     
#     #caculate deltaW numerically
#     deltaW = get_tone(Q, fs)
#     #print('Tone of Qf is {:e}'.format(deltaW))
#     
#     #mix to quadrature
#     QI, QQ = downmix(Q, fs, deltaW)
#     II, IQ = downmix(I, fs, deltaW)
#     
#     #convert data into real and complex components
#     cosAve = 2*(II-QQ)
#     sinAve = 2*(-QI-IQ)
#     
#     offset = np.arctan2(sinAve, cosAve)
#     #deltaW = get_tone(cosAve, fs)
#     #print('Tone of cosAve is {:e}'.format(deltaW))
#     
#     #phase rotate and calculate amps
#     cos2 = downmix(cosAve, fs, deltaW, single="Q")
#     sin2 = downmix(sinAve, fs, deltaW, single="I")
#     
#     amps = cos2 + sin2
#     
#     hyp = np.hypot(cosAve, sinAve)
#     
#     return amps, hyp, QQ, QI, IQ, II, deltaW
# =============================================================================
    
def DDC_ADC(signal, fs, ms, decim, int_time, output="real", diagnostics=False, finalFilt=False):
    
    """
    Rewritten, simplified DDC algorithm. Designed to closely approximate ADC EVM
    DDC settings.
    
    Inputs:
    signal: signal to be DDCed
    fs: sampling freqency
    ms: mixer frequency
    decim: how much to downsample the signal
    output: if "real" return real I and Q. if "complex", return complex data
    avgArray: if True, return averaged amplitude data over int_time length windows
    diagnostics: if true plot FFTs of signal at various stages of downconversion
    
    Outputs:
    data: if real, [I, Q]. if complex return a single data stream
    phase: phase reconstructed amplitude signal
    avgPhase: phase reconstructed averaged amplitudes
    """
    
    if diagnostics:
        FFT_plot([signal], fs, int_time, labels=["raw signal before DDC"])
        
    comp=signal
    
    #decimation factor for first decimation step
    decim1 = decim/2
    
    #decimation filter
    #passdB = .3 #maximum signal loss in pass band
    ntaps=128
    stopdB = 90 #minimum signal reduction in stop band
    passFrac = .8 #fraction of downsampled nyquist zone that is is in pass band
    dnyq = fs / (2*decim1) #width of Nyquist zone after downsampling
    #cutoffs = [ms - passFrac * dnyq, ms + passFrac*dnyq] #boundary of pass band
    passf = 1497e6 - 1480e6
    cutoffs = [passf, passf+6e6]
    transWidth = (1-passFrac) * dnyq #width of transition band
    
    a1 = 1
    b1 = firwin(ntaps, cutoffs, width=transWidth, window=('chebwin', stopdB), pass_zero=False, nyq=fs/2)
    if diagnostics:
        w, h = freqz(b1, a1, worN=3000000)
        compf = lfilter(b1, a1, comp)
        print(compf[0:5])
        fig, ax = FFT_plot([comp.real, compf.real], fs, int_time, labels=["I signal after mixing", "I signal after filtering"], returnFigAx=True)
        ax.plot(0.5*fs*w/(np.pi), 10*np.log10(np.abs(h)**2))
        ax.axvline(x=ms - passFrac * dnyq, label="Pass band cutoff", color='g')
        ax.axvline(x=ms + passFrac * dnyq, color='g')
        ax.axvline(x=ms - passFrac * dnyq - transWidth, label="Stop band cutoff", color='r')
        ax.axvline(x=ms + passFrac * dnyq + transWidth, color='r')
        ax.axvline(x=ms, label="mixer frequency", color="m")
        ax.axhline(y=-stopdB, label="Ideal supression in stop band", color='y')
        ax.set_title("Simulated LP filter response")
        ax.legend()
        fig.show()
        
        comp=compf
    else:
        comp = lfilter(b1, a1, comp)
        
    #downmix signal
    comp = downmix(signal, fs, ms, single="complex")
    
    if diagnostics:
        FFT_plot([comp.real], fs, int_time, labels=["I signal after mixing"])
    
    comp = decimate(comp, decim1)
    if diagnostics:
        FFT_plot([comp.real], fs/decim1, int_time, labels=["I signal after decimation"])
    
    #define filter parameters for final filter
    nyq=fs/2
    a=1
    b = firwin(ntaps-1, 2e7, width=1e7, nyq=nyq/decim1, pass_zero=False)
    comp=lfilter(b,a,comp)
    
    if diagnostics:
        w, h = freqz(b, a, worN=3000000)
        fig, ax = FFT_plot([comp.real], fs/decim1, int_time, labels=["I signal after LPF2"], returnFigAx=True)
        ax.plot(0.5*fs*w/(np.pi*decim1), 10*np.log10(np.abs(h)**2), label="LPF2 frequency response")
        ax.legend()
        fig.show()
    
    if output=="real":
        I = comp.real
        Q = 1j*comp.imag
        for i in range(0,4):
            I[i::4] = I[i::4]*(1j)**i
            Q[i::4] = Q[i::4]*(1j)**i
        I = I.real
        Q = Q.real
        if finalFilt==True:
            a2 = 1
            b2 = firwin(ntaps-1, 1e8, 1e8-7e7, nyq=nyq/decim1, pass_zero=False)
            I = lfilter(b2, a2, I)
            Q = lfilter(b2, a2, Q)
        if diagnostics:
            FFT_plot([I], fs/decim1, int_time, labels=["I signal after real conversion"])
        return I, Q, amplitude_avgs(comp, fs/decim, int_time)
        
    if output=="complex":
        #output DDCed data & averaged resolutions
        comp = decimate(comp, 2)
        if diagnostics:
            FFT_plot([comp.real], fs/decim1, int_time, labels=["after second decimation"])
        I = comp.real
        Q = comp.imag
        return I, Q, amplitude_avgs(comp, fs/decim, int_time)
        
def amplitude_avgs(signal, fs, int_time):
    """
    Calculates bin averaged amplitudes for a complex signal with hypotenuse reconstruction
    
    Inputs:
    signal: complex valued signal
    fs: sampling frequency
    int_time: window length to average over
    
    Outputs:
    avgs: average amplitude over each window
    """
    
    #figure out number of samples per window
    nsamps = int(fs * int_time)
    
    #determine number of windows to average over
    nwindows = int(np.floor(len(signal)/nsamps))
    
    #empty array to populate with averages
    avgs=[]
    
    #calculate averages
    for i in range(0, nwindows):
        window = signal[i:i+nsamps]
        hyps = np.hypot(window.real, window.imag)
        avg = np.mean(hyps)
        avgs.append(avg)
    
    return np.array(avgs)

def DDC_resolution(ampListA, ampListB, hist=True, scatter=True, title="", xlabel='Channel A asymmetry (PPM)', ylabel='Channel B asymmetry (PPM)'):
    """
    Calculates the amplitude asymmetry resolution from a list of bin averaged amplitudees.
    
    Inputs:
    ampListA: list of successive bin averaged amplitudes for channel A
    ampListB: list of successive bin averaged amplitudes for channel B
    hist: tells DDC_resolution whether or not to plot a resolution histogram
    scatter: tells DDC_resolution whether or not to plot a resolution scatter plot
    title: title for any plots produced
    
    Outputs:
    sigma: standard deviation in amplitude asymmetry difference
    mu: average amplitude asymmetry difference
    diffs: list of successive amplitude asymmetry diffs
    """
    
    #lists to store Channel A & B amplitude pairs
    Apairs = []
    Bpairs = []

    #populate Apairs & Bpairs with bin averaged RMS amplitudes
    for i in range(0, len(ampListA)-1, 2) :

        #start with A
        a1 = ampListA[i]
        a2 = ampListA[i+1]
        Apairs.append([a1, a2])

        #same process for B
        b1 = ampListB[i]
        b2 = ampListB[i+1]
        Bpairs.append([b1, b2])

    #convert lists to arrays to make calculations easier
    Apairs = np.array(Apairs)
    Bpairs = np.array(Bpairs)

    #calculate relative differences within channels
    Adiffs = (Apairs[:,0] - Apairs[:,1])/(Apairs[:,0] + Apairs[:,1])
    Bdiffs = (Bpairs[:,0] - Bpairs[:,1])/(Bpairs[:,0] + Bpairs[:,1])
    
    #calculate differences between channels
    diffs = Adiffs - Bdiffs

    #stat params
    sigma = np.std(diffs)
    mu = np.mean(diffs)
    
    #histogram plot the difference in amplitude asymmetry differences between channels
    if hist :
        f1, a1 = plt.subplots(1,1)
        y,x,_ = a1.hist(diffs*1e6, range = (-int(3*sigma*1e6), int(3*sigma*1e6)), bins=20)
        a1.set_xlabel('Amplitude Asymmetry Difference (PPM)')
        a1.text(-int(2.5*sigma*1e6), y.max(), s='Sigma={0:.3f} PPM'.format(sigma*1e6))
        a1.set_title(title)
        f1.show()
    
    #scatter plot
    if scatter :
        f2, a2 = plt.subplots(1,1)
        a2.scatter(Adiffs*1e6, Bdiffs*1e6, marker='.')
        a2.set_xlabel('Channel A asymmetry (PPM)')
        a2.set_ylabel('Channel B asymmetry (PPM)')
        a2.set_title(title)
        f1.show()

    return sigma, mu, diffs

def filter_finder(signal, fs, passdB, stopdB):
    """
    Recieves a signal with two frequency peaks and finds filter parameters that will
    filter out the higher frequency peak.
    
    Inputs:
    signal: the signal to find filter parameters for
    fs: sampling frequency of the signal
    passdB: maximum loss for the pass band
    stopdB: minimum loss for the stop band
    
    Outputs:
    IIRparams: parameters for an IIR filter, formated as IIRparams = [nyquist freq, pass freq, stop frequency, passdb, stopdb]
    """
    
    #frequency spectrum
    freqs, spectrum = periodogram(signal, fs=fs)
    
    #find the indices of peaks
    peaks, info = find_peaks(spectrum, height = .0000001)
    
    #find the peak frequencies
    peakFreqs = freqs[peaks[0]], freqs[peaks[1]]
    
    params = [np.max(freqs), peakFreqs[0], peakFreqs[1], passdB, stopdB]
    
    return params

def get_tone(signal, fs, res=False):
    """
    Returns the value of the largest frequency peak in a signal
    
    Inputs:
    signal: signal to be analyzed
    fs: sampling frequency of signal
    res: if True, return frequency resolution
    
    Outputs:
    tone: frequency the highest amplitude in the signal
    """
    
    #caculate the max frequency
    freqs, spectrum = periodogram(signal, fs=fs)
    peakInd = np.argmax(spectrum)
    tone = freqs[peakInd]
    
    if res:
        return tone, freqs[1]-freqs[0]
    else:
        return tone

def binToTxt(binPath, txtPath, bits, decimate) :
    """
    Takes a .bin file and converts it into MATLAB readable .txt file
    
    Inputs:
    binPath: location on disk of .bin file
    txtPath: location to save .txt file to (don't include .txt in this path)
    bits: resolution of incoming data (12 or 14 bits)
    decimate: if true, the .txt file that is saved will be 10x shorter than the .bin file (speeds up process to save data)
    """
    
    nchannels = 2 #break data into chA and chB. Don't further deinterleave into I & Q (for now - may change this later)
    midpoint = 2**(bits-1)
    data = np.fromfile(binPath, dtype=np.int16, count=-1, sep="")
    data = data-midpoint
    out=[]
    for i in range(nchannels) :
        out.append(data[i::nchannels])
    
    if decimate:
        out[0] = out[0][0:int(len(out[0])/10)]
        out[1] = out[1][0:int(len(out[1])/10)]
        
    out = 8*np.array(out).astype(int)
    
    np.savetxt(txtPath+'_chA.txt', out[0], fmt = '%i', delimiter='\n')
    np.savetxt(txtPath+'_chB.txt', out[1], fmt = '%i', delimiter='\n')
    
def signal_gen(fs, fo, amp, phi, offset, bits, noise, time, int_time, IQ) :
    """
    Returns a single frequency signal for testing DSP algorithms.
    
    Inputs:
        fs: sampling frequency
        fo: signal frequency
        amp: mid to peak amplitude of signal in units of volts (typical ADC value is .4-.6V)
        phi: intial phase of the signal
        offset: phase offset between channels
        bits: number of bits for output signal
        noise: magnitude of random noise in signal, in dBm (magnitude of noise is relative to amplitude of signal)
        time: length of the signal in time
        int_time: length of integration window
        IQ: boolean for whether or not to return I/Q components
        
    Outputs:
        A, B: signal output in channels A & B. If quad=true the return format is AI, AQ, BI, BQ
    """
    
    bins = np.int(np.floor(fs * time)) #number of bins for the length of the signal
    ratio = fo / fs #cycles per bin
    count = np.linspace(0,bins-1,bins) #array to store bin number values
    fftNum = int(fs)
    fsr = 1.35 #full scale range of ADC input channels in volts
    
    if IQ: #if quad=true return quadrature components of signal
        
        #AI = np.floor(amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi) + 2**(bits-1)*scale*np.random.rand(bins))
        #AQ = np.floor(amp*2**(bits-1)*np.cos(2*np.pi*ratio*count + phi) + 2**(bits-1)*scale*np.random.rand(bins))
        #BI = np.floor(amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi + offset) + 2**(bits-1)*scale*np.random.rand(bins))
        #BQ = np.floor(amp*2**(bits-1)*np.cos(2*np.pi*ratio*count + phi + offset) + 2**(bits-1)*scale*np.random.rand(bins))
        
        #generate the raw signal
        AIsig = amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi) / fsr
        AQsig = amp*2**(bits-1)*np.cos(2*np.pi*ratio*count + phi)/fsr
        BIsig = amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi + offset)/fsr
        BQsig = amp*2**(bits-1)*np.cos(2*np.pi*ratio*count + phi + offset)/fsr
        
        #calculate white noise parameters
        sigWatts = np.mean(AIsig**2)
        sigdB = 10 * np.log10(sigWatts)
        noisedB = sigdB + noise
        noiseWatts = 10**(noisedB/10)
        
        #scale noise amplitude to account for power density
        adj = int_time * fs / 2
        
        #generate white noise
        AInoiseAmp = np.random.normal(0, np.sqrt(adj*noiseWatts), size=bins)
        AQnoiseAmp = np.random.normal(0, np.sqrt(adj*noiseWatts), size=bins)
        BInoiseAmp = np.random.normal(0, np.sqrt(adj*noiseWatts), size=bins)
        BQnoiseAmp = np.random.normal(0, np.sqrt(adj*noiseWatts), size=bins)
        
        #comine noise and signal
        AI = np.floor(np.clip(AInoiseAmp+AIsig, -2**(bits-1), 2**(bits-1)))
        AQ = np.floor(np.clip(AQnoiseAmp+AQsig, -2**(bits-1), 2**(bits-1)))
        BI = np.floor(np.clip(BInoiseAmp+BIsig, -2**(bits-1), 2**(bits-1)))
        BQ = np.floor(np.clip(BQnoiseAmp+BQsig, -2**(bits-1), 2**(bits-1)))
        
        return AI, AQ, BI, BQ
    
    else : #just return the real valued signals
        
        #A = np.floor(amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi) + amp*2**(bits-1)*np.random.normal(0, 10**(noise/20), size=bins))
        #B = np.floor(amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi + offset) + amp*2**(bits-1)*np.random.normal(0, 10**(noise/20), size=bins))
        
        #generate the raw signal
        Asig = amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi)/fsr
        Bsig = amp*2**(bits-1)*np.sin(2*np.pi*ratio*count + phi + offset)/fsr
        
        #calculate white noise parameters
        sigWatts = np.mean(Asig**2)
        sigdB = 10 * np.log10(sigWatts)
        noisedB = sigdB + noise
        noiseWatts = 10**(noisedB/10)
        
        #scale noise amplitude to account for power density
        adj = int_time * fs / 2
        
        #generate white noise
        AnoiseAmp = np.random.normal(0, np.sqrt(adj*noiseWatts), size=bins)
        BnoiseAmp = np.random.normal(0, np.sqrt(adj*noiseWatts), size=bins)
        
        #comine noise and signal
        A = np.floor(np.clip(AnoiseAmp+Asig, -2**(bits-1), 2**(bits-1)))
        B = np.floor(np.clip(BnoiseAmp+Bsig, -2**(bits-1), 2**(bits-1)))
        
        return A, B

def bin_avg(signal, fs, int_time):
    """
    Similar to amplitude_avgs, but returns the mean value across a series of bins for a given signal.
    
    Inputs:
    signal: real valued signal
    fs: sampling frequency
    int_time: window length to average over
    
    Outputs:
    avgs: average value over each window
    """
    
    #figure out number of samples per window
    nsamps = int(fs * int_time)
    
    #determine number of windows to average over
    nwindows = int(np.floor(len(signal)/nsamps))
    
    #empty array to populate with averages
    avgs=[]
    
    #calculate averages
    for i in range(0, nwindows):
        window = signal[i:i+nsamps]
        hyps = np.hypot(window.real, window.imag)
        avg = np.mean(hyps)
        avgs.append(avg)
    
    return np.array(avgs)

def hypot_recon(AI, AQ, BI, BQ, fs, int_time, raw):
    """
    Hypotenuse amplitude reconstruction
    Inputs:
        AI, AQ, BI, BQ: in phase and quadrature components of channel A & B signals
        fs: sampling frequency
        raw: if True, returns raw amplitude signal, otherwise binned amplitude averages
    Outputs:
        ampA, ampB: reconstructed amplitude for the signals from A&B
    """
    
    #pointwise ampmlitude calculation
    rawAmpA = np.hypot(AI, AQ)
    rawAmpB = np.hypot(BI, BQ)
    
    if raw:
        return rawAmpA, rawAmpB
    else:
        return bin_avg(rawAmpA, fs, int_time), bin_avg(rawAmpB, fs, int_time)

def phase_rotation(AI, AQ, BI, BQ, fs, fo, int_time, raw) :
    """
    Phase rotation amplitude reconstruction
    Inputs:
        AI, AQ, BI, BQ: in phase and quadrature components of channel A & B signals
        fs: sampling frequency
        fo: signal frequency
        int_time: window length
        raw: if True, returns raw amplitude signal, otherwise binned amplitude averages
    Outputs:
        ampA, ampB: reconstructed amplitude for the signals from A&B
    """
    
    #mix the signals down into II, QQ, QI, and IQ components
    ratio = fo/fs
    count = np.linspace(1, len(AI), len(AI))
    I = np.sin(2*np.pi*ratio*count) #in phase mixer
    Q = np.cos(2*np.pi*ratio*count) #quadrature mixer
    AII = AI*I
    AIQ = AI*Q
    AQI = AQ*I
    AQQ = AQ*Q
    BII = BI*I
    BIQ = BI*Q
    BQI = BQ*I
    BQQ = BQ*Q
    
    #cancel out high frequency components, hold onto zero-frequency components
    AxPhase = AII + AQQ
    AyPhase = AIQ - AQI
    BxPhase = BII + BQQ
    ByPhase = BIQ - BQI
    
    #calculate pointwise amplitude signals
    rawAmpA = np.hypot(AxPhase, AyPhase)
    rawAmpB = np.hypot(BxPhase, ByPhase)
    
    if raw:
        return rawAmpA, rawAmpB
    else:
        #calculate and return bin averaged amplitudes
        return bin_avg(rawAmpA, fs, int_time), bin_avg(rawAmpB, fs, int_time)

def binned_phase_rotation(AI, AQ, BI, BQ, fs, int_time, raw) :
    """
    Phase rotation amplitude reconstruction but without foreknowledge of fo.
    Takes FFT on each bin and uses that for amplitude reconstruction
    
    Inputs:
        AI, AQ, BI, BQ: in phase and quadrature components of channel A & B signals
        fs: sampling frequency
    Outputs:
        ampA, ampB: binned reconstructed amplitude for the signals from A&B
    """
    
    #figure out number of samples per window
    nsamps = int(fs * int_time)
    
    #determine number of windows to average over
    nwindows = int(np.floor(len(AI)/nsamps))
    
    #empty arrays to populate with bin averages
    Aavgs=[]
    Bavgs=[]
    Araw=[]
    Braw=[]
    
    #calculate averages
    for i in range(0, nwindows):
        
        #isolate the given window
        windowAI = AI[i:i+nsamps]
        windowAQ = AQ[i:i+nsamps]
        windowBI = BI[i:i+nsamps]
        windowBQ = BQ[i:i+nsamps]
        
        #calculate local max frequency
        toneAI = get_tone(windowAI, fs)
        toneAQ = get_tone(windowAQ, fs)
        toneBI = get_tone(windowBI, fs)
        toneBQ = get_tone(windowBQ, fs)
        toneA = np.mean([toneAI, toneAQ])
        toneB = np.mean([toneBI, toneBQ])
        
        #mix the signals down into II, QQ, QI, and IQ components
        ratioA = toneA/fs
        ratioB = toneB/fs
        count = np.linspace(0, nsamps-1, nsamps)
        
        #create mixer signals
        IA = np.sin(2*np.pi*ratioA*count) #in phase mixer
        QA = np.cos(2*np.pi*ratioA*count) #quadrature mixer
        IB = np.sin(2*np.pi*ratioB*count)
        QB = np.cos(2*np.pi*ratioB*count)
        
        #mix
        AII = windowAI*IA
        AIQ = windowAI*QA
        AQI = windowAQ*IA
        AQQ = windowAQ*QA
        BII = windowBI*IB
        BIQ = windowBI*QB
        BQI = windowBQ*IB
        BQQ = windowBQ*QB
        
        #cancel out high frequency components, hold onto zero-frequency components
        AxPhase = AII + AQQ
        AyPhase = AIQ - AQI
        BxPhase = BII + BQQ
        ByPhase = BIQ - BQI
        
        #calculate pointwise amplitude signals
        rawAmpA = np.hypot(AxPhase, AyPhase)
        rawAmpB = np.hypot(BxPhase, ByPhase)
        
        #average over the bins and store
        Araw.append(rawAmpA)
        Braw.append(rawAmpB)
        Aavgs.append(np.mean(rawAmpA))
        Bavgs.append(np.mean(rawAmpB))
    
    #return bin amplitudes
    if raw:
        return Araw, Braw
    else:
        return Aavgs, Bavgs
    
# =============================================================================
# def fourier_amp(A, B, fs) :
#     """
#     Determines signal amplitude from magnitude of largest peak in frequency space
#     Inputs:
#         A, B: Signals from channels A and B
#     Outputs:
#         Aamp, Bamp: continuously reconstructed amplitude for the signals from A&B
#     """
#     
#     #find the max frequency
#     Afreqs, Aspectrum = periodogram(A, fs=fs)
#     ApeakInd = np.argmax(Aspectrum)
#     Atone = Afreqs[ApeakInd]
#     Aamp = np.abs(Aspectrum[ApeakInd])
#     
#     #find the max frequency
#     Bfreqs, Bspectrum = periodogram(B, fs=fs)
#     BpeakInd = np.argmax(Bspectrum)
#     Btone = Bfreqs[BpeakInd]
#     Bamp = np.abs(Bspectrum[BpeakInd])
#     
#     return Aamp, Bamp
# =============================================================================

def least_sq(sig1, sig2) :
    """
    Calculates and returns the least squares difference between two signals.
    Inputs:
        sig1/sig2 : signals to have their least squares difference returned
    Outputs:
        diff: the least squares diffrence
    """
    return np.sum(np.sqrt(np.abs(sig1-sig2)))

def amplitude_correlation(signal1, signal2, title='', xlabel='',ylabel='', hist=True, scatter=False):
    """
    Measures the correlation between a list of amplitudes from a signal and its corresponding "truth" values
    Inputs:
        signal1/signal2: list of amplitudes
        scatter: instruction for making a correlation scatterplot
        hist: instruction for making a histogram
        title, xlabel, ylabel: labels for scatter plot and histogram
    Outputs:
        sigma: standard deviation of differenceasymmetry between signals
        asym: binwwise asymmetry between the input signals
    """
    
    #make signals into arrays if they aren't already
    signal1 = np.array(signal1)
    signal2 = np.array(signal2)
    
    #calculate asymmetry and sigma
    asym = (signal1 - signal2) / (signal1 + signal2)
    sigma = np.std(asym)
    mu = np.mean(asym)
    
    if hist :
        
        sigLab = matplotlib.offsetbox.AnchoredText('Sigma={0:.3f} PPM'.format(sigma*1e6), loc=2)
        f1, a1 = plt.subplots(1,1)
        y,x,_ = a1.hist(asym*1e6, bins=15)
        a1.set_xlabel('Asymmetry (ppm)')
        a1.add_artist(sigLab)
        a1.set_title(title)
        f1.show()
        
    if scatter :
        sigLab = matplotlib.offsetbox.AnchoredText('Sigma={0:.3f} PPM'.format(sigma*1e6), loc=2)
        f2, a2 = plt.subplots(1,1)
        a2.scatter(signal1, signal2, marker='.')
        a2.plot(signal1, signal2, color='none')
        a2.relim()
        a2.autoscale_view()
        a2.set_xlabel(xlabel)
        a2.set_ylabel(ylabel)
        a2.set_title(title)
        a2.add_artist(sigLab)
        f2.show()
    
    return sigma, asym

def power_fitter(xdata, ydata, plot=True, title="", dataLabel = "") :
    """
    Creates a log-log plot, finds and returns the best fit slope of a power law line (y = a x^b since log-log)
    Inputs:
        xdata: x axis info
        ydata: y axis info
        
    Outputs:
        a, b - intercept and slope of best fit power function
    """
    
    #find the least squares, polynomial fit to the log of the data, rather than data itself
    #accounts for logarithmic growth without over fitting to large number data
    #and returns correct parameters for power law y = a x^b
    ab = np.polyfit(np.log(xdata),np.log(ydata), 1)
    a = np.exp(ab[1])
    b = ab[0]
    
    if plot:
        comp = matplotlib.offsetbox.AnchoredText('Best fit time complexity: O(n^{:.3f})\nBest fit power equation: y = {:.3e} x ^ {:.3f}'.format(b,a,b), loc=1)
        fig, ax = plt.subplots(1,1)
        ax.plot(xdata, ydata, label=dataLabel)
        ax.plot(xdata, a*xdata**b, label="Best fit power law line")
        ax.set_title(title)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.add_artist(comp)
        ax.legend(loc=2)
        fig.show()
    
    return a, b

def exp_fitter(xdata, ydata, plot=True, title="", dataLabel = "") :
    """
    Creates a log-log plot, finds and returns the best fit slope of an exponential
    Inputs:
        xdata: x axis info
        ydata: y axis info
        
    Outputs:
        a, b - intercept and slope of best fit exponential function y=ae^(bx)
    """
    
    ab = np.polyfit(xdata,np.log(ydata), 1)
    a = np.exp(ab[1])
    b = ab[0]
    
    if plot:
        comp = matplotlib.offsetbox.AnchoredText('Best fit exponential equation: y = {:.3e} e^({:.3e} x)'.format(a,b), loc=1)
        fig, ax = plt.subplots(1,1)
        ax.loglog(xdata, ydata, label=dataLabel)
        ax.loglog(xdata, a*np.e**(b*xdata), label="Best fit exponential")
        ax.set_title(title)
        ax.add_artist(comp)
        ax.legend(loc=2)
        fig.show()
    
    return a, b

def quasiPoly_fitter(xdata, ydata, plot=True, title="", dataLabel = "") :
    """
    Creates a log-log plot, finds and returns the best fit slope of an quasi-polynomial growth equation
    Inputs:
        xdata: x axis info
        ydata: y axis info
        
    Outputs:
        a, b - best fit parameters for function y=ax^(b log(x))
    """
    
    ab = np.polyfit(np.log(xdata)**2,np.log(ydata), 1)
    a = np.exp(ab[1])
    b = ab[0]
    
    if plot:
        comp = matplotlib.offsetbox.AnchoredText('Best fit quasi-polynomial time equation: y = {:.3e} x^({:.3e} log(x))'.format(a,b), loc=1)
        fig, ax = plt.subplots(1,1)
        ax.plot(xdata, ydata, label=dataLabel)
        ax.plot(xdata, a*xdata**(b*np.log(xdata)), label="Best fit quasi-polynomial line")
        ax.set_title(title)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.add_artist(comp)
        ax.legend(loc=2)
        fig.show()
    
    return a, b
    