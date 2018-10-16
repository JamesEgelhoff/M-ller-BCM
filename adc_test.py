import tools as t
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fileraw = "C:/Users/molle/Desktop/ADC32RF45/data/8_3_2018_1400MA_3000MS_12bit_200ms.bin"
filevoltage = "C:/Users/molle/Desktop/ADC32RF45/data/8_3_2018_1400MA_3000MS_12bit_200ms.npy"
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

t.read_binary(infile=fileraw, outfile=filevoltage, bits=12, fsr=1.35, raw=False)
A, B = t.open_binary(filevoltage)
# t.xcor_spectrum(A, B, fo, fs, int_time=int_time, n_window = ncalc, dual=True)
#resList, widthList = t.plot_res_vs_binwidth(A, B, 1e-6, 1.001e-3, 1001, log=True)
# sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, bits=bits, ncalc = ncalc, pulse_width=int_time, hist=False, scatter=True)

plt.plot(A[0:500])
plt.plot(B[0:500])
plt.show()

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
