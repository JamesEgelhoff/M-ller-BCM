import tools as t
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import iirdesign, freqz, lfilter

fileraw = "C:/Users/molle/Desktop/ADC32RF45/data/11_5_2018_1497MA_8dB_3000MS_12bit_200ms_DDCx4_1450NCO.bin"
filevoltage = "C:/Users/molle/Desktop/ADC32RF45/data/11_5_2018_1497MA_8dB_3000MS_12bit_200ms_DDCx4_1450NCO.npy"
# fileavg = '../ADC_DATA/8_3/1400_3000_100_avg_amplitude.npy'
fo = 1.497e9
foff = -38720 # 10 Hz accuracy (calculated with 100ms fft)
fs=.75e9
int_time = .5*1e-3
fsr= 1.35 # ADC32RF45 fullscale range (volts)
bits = 16
ncalc = 200
filter_ntaps = 1024
npt2n = 2**29 #= 536870912, this value is around 180ms around 3GHz sampling.
nchannels = 4
title = "Setup with new Mini-Circuits parts and x4 DDC"
ms = 1450e6 #mixer frequency for DDC calcs

DDC = True
if DDC:
    dtype = np.uint16
else:
    dtype = np.int16

#t.read_binary(infile=fileraw, outfile=filevoltage, bits=bits, fsr=1.35, send='mid', dtype=dtype, nchannels=nchannels)
#A, B, C, D = t.open_binary(filevoltage, nchannels=nchannels)
#t.xcor_spectrum(A,C, fo, fs, int_time=int_time, n_window = "auto", dual=True, dualTitle=title)
#resList, widthList = t.plot_res_vs_binwidth(A, B, 1e-6, 1.001e-3, 5001, log=True, title=title)
#sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, ncalc = "auto", pulse_width=int_time, hist=True, scatter=True, title=title)
#ChA_amp, ChB_amp = t.DDC_amplitudes(A, B, C, D, plot=True, title=title)


#code to test DDC_amplitudes function
#make a reference signal, mix with a sin and cos signal to produce I & Q. Filter and test amplitude reconstruction
sig_phi = np.random.rand()*2*np.pi # random phase for signal
sig_phi_off = .2*np.random.rand()*np.pi #small phase offset between two simulated channels
sig_f = 1497e6 #signal freq, Hz
ref_f = 1503e6 #signal due to reflection over nyquist boundary
mixer = 1450e6 #mixer freq
samples = int(6e7) # number of samples to work with
fs = 3e9 #sampling freqency
time = np.linspace(0, samples/fs, samples)

#make the signals for Channels A and B
A_signal =  .5*(np.sin(time*sig_f+sig_phi) + np.sin(time*ref_f+sig_phi))
A_I = A_signal*np.sin(time*mixer)
A_Q = A_signal*np.cos(time*mixer)

B_signal =  .5*(np.sin(time*sig_f+sig_phi+sig_phi_off)+np.sin(time*ref_f+sig_phi+sig_phi_off))
B_I = B_signal*np.sin(time*mixer)
B_Q = B_signal*np.cos(time*mixer)

#LP filter the downmixed signals w/ IIR
cutoff=2*np.abs(mixer - sig_f)
nyq = 0.5 * fs
normalized_pass = cutoff / nyq
normalized_stop = 100*cutoff / nyq
b0, a0 = iirdesign(normalized_pass, normalized_stop, .3, 60)
w, h = freqz(b0, a0, worN=3000000)
A_If = lfilter(b0, a0, A_I)
A_Qf = lfilter(b0, a0, A_Q)
B_If = lfilter(b0, a0, B_I)
B_Qf = lfilter(b0, a0, B_Q)

#useful diagnostic info
fig, ax = plt.subplots(1,1)
ax.plot(A_signal[0:1000], label="raw, unmixed simulated data")
ax.plot(A_If[0:1000], label = "simulated downmixed and filtered data, ChA I stream")
ax.plot(A_Qf[0:1000], label = "simulated downmixed and filtered data, ChA Q stream")
ax.legend(loc=3)
plt.show()

fig2, ax2 = plt.subplots(1,1)
ax2.plot(0.5*fs*w/np.pi, 10*np.log10(np.abs(h)))
ax2.set_title('Filter response (dBm)')
ax2.set_xlabel('Frequency (Hz)')
ax2.set_ylabel('Gain (dBm)')
ax2.set_xscale('log')
fig2.show()

t.xcor_spectrum(A_If,B_If, fo, fs, int_time=int_time, n_window = "auto", dual=True, dualTitle="Simulated data spectrum, downmixed and filtered, including reflected signal")

#test DDC_amplitudes
ChA_amp, ChB_amp = t.DDC_amplitudes(A_If, A_Qf, B_If, B_Qf, plot=True, title="Hypotenuse amplitude reconstruction on simulated downmixed and filtered data, including reflected signal")


# =============================================================================
# plt.plot(A[0:700])
# plt.plot(B[0:700])
# plt.title(title)
# plt.show()
# =============================================================================
    
# =============================================================================
# tone_f, foff = t.get_freq(data=A, fs=fs, fo_nominal=fo, int_time = int_time, plot=False)
# mixoff=.1*np.pi*foff
# 
# a, avg, avg2, raw_avgs = t.ddc(ChA=A, ChB=B, fo=fo, foff=foff, mixoff=mixoff, lpf_fc=np.abs(mixoff), lpf_ntaps=filter_ntaps, 
#                      fs = fs, bits = bits, int_time=int_time, ncalc=ncalc, calc_off=filter_ntaps+10,
#                      phase_time=int_time, nch=2, plot_en=True, plot_len=3e6, plot_win=0,
#                      plot_Fourier=True, filt="cheby1", suppress=False, amp_res=True)
# =============================================================================

input('press enter to finish ADC script.')
