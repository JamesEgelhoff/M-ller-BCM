# math
import numpy as np
from scipy.fftpack import fft, rfft
from scipy.integrate import trapz
from scipy.signal import medfilt, butter, bessel, firwin, lfilter, freqz, decimate, periodogram, welch, iirdesign, get_window, firls
from scipy.stats import norm

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

def read_binary(infile, outfile, bits=12, fsr=1.35, raw=None) :
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
    if raw :
        out = np.fromfile(infile, dtype=np.int16, count=-1, sep="")
    else :
        data = (np.fromfile(infile, dtype=np.int16, count=-1, sep="")-(midpoint))*fsr/(maxcode)
        out = np.stack((data[0::2], data[1::2]))
    np.save(outfile, out)

def open_binary(infile) :
    data = np.load(infile)
    return data[0], data[1]


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

    input('press enter to finish frequency script.')
    
    return tone_f, f_off

def xcor_spectrum(ChA, ChB, fo, fs, nbits=12, int_time=1e-3, n_window=99, plot=False, dual=False, dualTitle="") :
    '''calculate the fourier spectrum of the adc's digitized signal.
    convert spectrum to units of dBc/Hz (decibels normalized to the carrier power in 1Hz bandwidth).
    calculate the jitter contribution in some specified bandwidth relative to the carrier.
    input:
    file- string giving binary file containing data
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

    for i in range(n_window):
        print(i)
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
        plt.step(fbins, 10*np.log10(Saa/np.max(Saa)), label='CH A Noise Spectral Density')
        plt.step(fbins, 10*np.log10(Sbb/np.max(Sbb)), label='CH B Noise Spectral Density')
        plt.xscale('log')
        plt.xlabel('Frequency (Hertz)')
        plt.ylabel(r'$\frac{dBc}{Hz}$', rotation=0, fontsize=16)
        plt.legend(loc=2)
        plt.title(dualTitle)
        plt.show()

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
    plt.plot(0.5*fs*w/np.pi, np.abs(h), 'b')
    plt.plot(np.abs(mixoff), 0.5*np.sqrt(2), 'ko')
    plt.axvline(cutoff, color='k')
    if filt=="IIR":
        plt.axvline(normalized_stop*nyq, color='k')
    #plt.xlim(0, 0.5*fs)
    plt.title("Lowpass Filter Frequency Response")
    plt.xlabel('Frequency [Hz]')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.show()

    # arrays to store average amplitudes
    avg = np.zeros((nch,ncalc))
    avg2 = np.zeros((nch,ncalc))
    
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

    return a, avg, avg2

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
    
def amplitude_asym_hist(ChA, ChB, fs=3e9, bits=12, ncalc="auto", pulse_width=5e-4, hist=True, scatter=False, title=""):
    
    
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
        
        print(i)

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
        y,x,_ = plt.hist(diffs*1e6, range = (-int(3*sigma*1e6), int(3*sigma*1e6)), bins=20)
        plt.xlabel('Amplitude Asymmetry Difference (PPM)')
        plt.text(-int(2.5*sigma*1e6), y.max(), s='Sigma={0:.3f} PPM'.format(sigma*1e6))
        plt.title(title)
        plt.show()
    
    #scatter plot
    if scatter :
        plt.scatter(Adiffs*1e6, Bdiffs*1e6)
        plt.xlabel('Channel A asymmetry (PPM)')
        plt.ylabel('Channel B asymmetry (PPM)')
        plt.title(title)
        plt.show()

    return sigma, mu, diffs

def plot_res_vs_binwidth(ChA, ChB, minWidth, maxWidth, steps, title, ncalc = 30, fs=3e9, bits=12, log=False) :
    
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
        print('Calculating resolution for {0:f} millisecond pulse windows'.format(1e3*widthList[i]))
        res, mu, diffs = amplitude_asym_hist(ChA, ChB, fs=fs, bits=bits, ncalc=ncalc, pulse_width=widthList[i], hist=False, scatter=False)
        resList[i] = res
    
    if log==False :
        plt.plot(widthList*1e3, resList*1e6, '.')
        plt.xlabel('Pulse width (ms)')
        plt.ylabel('Resolution (PPM)')
        plt.title(title)
        plt.show()
    else:
        plt.plot(widthList*1e3, resList*1e6, '.')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Pulse width (ms)')
        plt.ylabel('Resolution (PPM)')
        plt.title(title)
        plt.show()
    
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