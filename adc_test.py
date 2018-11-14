import tools as t
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import iirdesign, freqz, lfilter

fileraw = "C:/Users/molle/Desktop/ADC32RF45/data/11_13_2018_1330MA_10dB_3000MS_12bit_200ms.bin"
filevoltage = "C:/Users/molle/Desktop/ADC32RF45/data/11_13_2018_1330MA_10dB_3000MS_12bit_200ms.npy"
# fileavg = '../ADC_DATA/8_3/1400_3000_100_avg_amplitude.npy'
fo = 1.497e9
#fo=1330e6
foff = -38720 # 10 Hz accuracy (calculated with 100ms fft)
fs=3e9
int_time = .5*1e-3
fsr= 1.35 # ADC32RF45 fullscale range (volts)
ncalc = 200
filter_ntaps = 1024
npt2n = 2**29 #= 536870912, this value is around 180ms around 3GHz sampling.
title = ""
ms = 1450e6 #mixer frequency for DDC calcs

#DDC settings.
DDC = False
bits=12
band = "single"
decim=8
if DDC:
    bits=16
    dtype = np.uint16
    fs=fs/decim
    nchannels=4
else:
    dtype = np.int16
    nchannels=2

#open files
#t.read_binary(infile=fileraw, outfile=filevoltage, bits=bits, fsr=1.35, send='mid', dtype=dtype, nchannels=nchannels)
if DDC==False:
    A, B = t.open_binary(filevoltage, nchannels=nchannels)
else:
    AQ, AI, BI, BQ = t.open_binary(filevoltage, nchannels=nchannels)

#t.xcor_spectrum(A, C, fo, fs, int_time=int_time, n_window = "auto", dual=True, dualTitle=title)
#resList, widthList = t.plot_res_vs_binwidth(A, B, 1e-6, 1.001e-3, 5001, log=True, title=title)
# sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, ncalc = "auto", pulse_width=int_time,
#                                          hist=True, scatter=True, 
#                                          title="Resolution w/ 100ms 1497MHz signal, no DDC")
#ChA_amp, ChB_amp = t.DDC_amplitudes(A, B, C, D, plot=True, title=title)

A = A[0:int(1.5e8)]
B = B[0:int(1.5e8)]

sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, ncalc = "auto", pulse_width=int_time,
                                         hist=True, scatter=True, 
                                         title="Resolution w/ 100ms 1330MHz signal, no DDC")

# =============================================================================
# fig, ax = plt.subplots(1,1)
# ax.plot(A[0:800], label="Channel A")
# ax.plot(hypA[0:800], label="hypotenuse reconstructed A amplitude after DDC")
# ax.legend(loc=3)
# ax.set_xlabel("sample number")
# ax.set_ylabel("amplitude in volts")
# fig.show()
# =============================================================================


input('press enter to finish ADC script.')
