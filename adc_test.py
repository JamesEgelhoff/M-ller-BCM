import tools as t
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import iirdesign, freqz, lfilter, periodogram, find_peaks

fileraw = "E:/more ADC data/11_27_2018_1330MA_13dB_3000MS_200ms_DDCx16_1300NCO_real.bin"
filevoltage = "E:/more ADC data/11_27_2018_1330MA_13dB_3000MS_200ms_DDCx16_1300NCO_real.npy"
# fileavg = '../ADC_DATA/8_3/1400_3000_100_avg_amplitude.npy'

fs=3e9
fo = 1.497e9
foff = -38720 # 10 Hz accuracy (calculated with 100ms fft)
int_time = .5*1e-3
fsr= 1.35 # ADC32RF45 fullscale range (volts)
ncalc = 200
filter_ntaps = 32
npt2n = 2**29 #= 536870912, this value is around 180ms around 3GHz sampling.
title = ""
ms = 1300e6 #mixer frequency for DDC calcs

#DDC settings.
DDC = True
bits=12
band = "single"
decim=16
if DDC:
    bits=16
    dtype = np.uint16
    fs=fs/decim
    nchannels=4
else:
    dtype = np.int16
    nchannels=2

# =============================================================================
# #open files
# t.read_binary(infile=fileraw, outfile=filevoltage, bits=bits, fsr=1.35, send='mid', dtype=dtype, nchannels=nchannels)
# if DDC==False:
#     A, B = t.open_binary(filevoltage, nchannels=nchannels)
# else:
#     AI, AQ, BI, BQ = t.open_binary(filevoltage, nchannels=nchannels)
# =============================================================================

#t.xcor_spectrum(A, C, fo, fs, int_time=int_time, n_window = "auto", dual=True, dualTitle=title)
#resList, widthList = t.plot_res_vs_binwidth(A, B, 1e-6, 1.001e-3, 5001, log=True, title=title)
# sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, ncalc = "auto", pulse_width=int_time,
#                                          hist=True, scatter=True, 
#                                          title="Resolution w/ 100ms 1497MHz signal, no DDC")
#ChA_amp, ChB_amp = t.DDC_amplitudes(A, B, C, D, plot=True, title=title)
    
A = A[0:int(6e7)]
fo=t.get_tone(A, fs)
#B=None
fs=3e9
ntaps=32
#I, Q, avgs = t.DDC_ADC(A, fs, ms, decim, int_time, output="real", diagnostics=False)

title = "I and Q for ADC32RF45 DDC, 1330MHz signal, 1300MHz mixer"
fig, ax = plt.subplots(1,1)
#ax.plot(AI[0:500], label="I data stream from hardware DDC")
#ax.plot(AQ[0:500], label="Q data stream from hardware DDC")
ax.plot(I[ntaps:500+ntaps], label="I data stream from software DDC")
ax.plot(Q[ntaps:500+ntaps], label="Q data stream from software DDC")
ax.set_title(title)
ax.set_ylabel("Voltage")
ax.legend()
fig.show()

t.FFT_plot([I, AI], fs/decim, int_time, labels=["simulated in phase component", "hardware in phase component"])

input('press enter to finish ADC script.')
