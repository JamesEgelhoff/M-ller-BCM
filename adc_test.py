import tools as t
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileraw = "C:/Users/molle/Desktop/ADC32RF45/data/10_17_2018_1497MA_15dB_3000MS_12bit_200ms_setup6.bin"
filevoltage = "C:/Users/molle/Desktop/ADC32RF45/data/10_17_2018_1497MA_15dB_3000MS_12bit_200ms_setup6.npy"
=======
fileraw = "C:/Users/James/Desktop/Moller/ADC data/10_17_2018_1497MA_15dB_3000MS_12bit_200ms_setup6.bin"
filevoltage = "C:/Users/James/Desktop/Moller/ADC data/10_17_2018_1497MA_15dB_3000MS_12bit_200ms_setup6.npy"
>>>>>>> 441095fc4a5e31eccae6e6c090d28bd451fa37a7
# fileavg = '../ADC_DATA/8_3/1400_3000_100_avg_amplitude.npy'
fo = 1.497e9
foff = -38720 # 10 Hz accuracy (calculated with 100ms fft)
fs=3e9
int_time = .5*1e-3
fsr= 1.35 # ADC32RF45 fullscale range (volts)
bits = 12
ncalc = 400
filter_ntaps = 100
npt2n = 2**29 #= 536870912, this value is around 180ms around 3GHz sampling.
title="Signal Power vs Resolution in Setup 6 with 20ms data captures"
adcZ = 65 #ohms

# =============================================================================
# #little script to process data taken at various signal amplitudes
# fileList = []
# for i in range(-5, 16):
#     if i < 0:
#         binString = "C:/Users/molle/Desktop/ADC32RF45/data/10_18_2018_1497MA_m{:}dB_3000MS_12bit_20ms.bin".format(-1*i)
#         npyString = "C:/Users/molle/Desktop/ADC32RF45/data/10_18_2018_1497MA_m{:}dB_3000MS_12bit_20ms.npy".format(-1*i)
#     else:
#         binString = "C:/Users/molle/Desktop/ADC32RF45/data/10_18_2018_1497MA_{:}dB_3000MS_12bit_20ms.bin".format(i)
#         npyString = "C:/Users/molle/Desktop/ADC32RF45/data/10_18_2018_1497MA_{:}dB_3000MS_12bit_20ms.npy".format(i)
#         
#     t.read_binary(infile=binString, outfile=npyString, bits=bits, fsr=fsr, raw=False)
#     fileList.append(npyString)
#     print(i)
# =============================================================================

powerResList = t.power_scatter(fileList, adcZ, title=title)

#t.read_binary(infile=fileraw, outfile=filevoltage, bits=12, fsr=1.35, raw=False)
#A, B = t.open_binary(filevoltage)
#t.xcor_spectrum(A, B, fo, fs, int_time=int_time, n_window = ncalc, dual=True, dualTitle = title)
#resList, widthList = t.plot_res_vs_binwidth(A, B, 1e-6, 1.001e-3, 1001, title, log=True)
#sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, bits=bits, ncalc = "auto", pulse_width=int_time, hist=True, scatter=False, title=title)

# =============================================================================
# plt.plot(A[0:500], label="Channel A")
# plt.plot(B[0:500], label="Channel B")
# plt.legend(loc=2)
# plt.title(title)
# plt.show()
# =============================================================================

# freqAmp = np.array([MA1497_m3dB*1e6, MA1497_0dB*1e6, MA1497_3dB*1e6, MA1497_6dB*1e6, MA1497_9dB*1e6])

# =============================================================================
# fig = plt.figure()
# ax1 = fig.add_subplot(111,projection='3d')
# 
# lx = len(freqAmp[0])
# ly = len(freqAmp[:,0])
# xpos = np.arange(1197, 1197+100*lx, 100)
# ypos = np.arange(3, 3+3*ly, 3)
# xpos, ypos = np.meshgrid(xpos, ypos)
# 
# xpos = xpos.flatten()
# ypos = ypos.flatten()
# zpos = np.zeros(lx*ly)
# 
# dx = 30*np.ones_like(zpos)
# dy = .7*np.ones_like(zpos)
# dz = freqAmp.flatten()
# 
# ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, shade=True)
# ax1.set_xlabel("Frequency(MHz)")
# ax1.set_ylabel("Carrier Amplitude (dB)")
# ax1.set_zlabel("Relative uncertainty")
# ax1.set_title('Signal resolution vs. Amplitude and Frequency')
# ax1.ticklabel_format(style='sci', axis='z')
# =============================================================================

# =============================================================================
# err = freqAmp/np.sqrt(400)
# plt.errorbar([-3, 0, 3, 6, 9], freqAmp, fmt='o', yerr = err)
# plt.xlabel('Signal Amplitude (dB)')
# plt.ylabel('Resolution (PPM)')
# plt.show()
# =============================================================================

# foff = t.get_freq(data=A, fs=fs, fo_nominal=fo, int_time = 10e-3, plot=False)

#a, avg, avg2 = t.ddc(ChA=A, ChB=B, fo=fo+foff, lpf_fc=1e3, lpf_ntaps=filter_ntaps, fs = fs, bits = bits, int_time=int_time,
#                     ncalc=ncalc, calc_off=filter_ntaps+10, phase_time=int_time, nch=2, plot_en=False, plot_len=3e6, plot_win=0)

#print('amplitude distribution =', np.std(avg[0]/avg[1]), 'sigma (phase reconstruction)\n')
#print('amplitude distribution =', np.std(avg2[0]/avg2[1]), 'sigma (hypotenuse reconstruction)')



# s.res_plot(avg2[0], avg2[1])

# np.save(fileavg, avg2)
# data = np.load(fileavg)

input('press enter to finish ADC script.')
