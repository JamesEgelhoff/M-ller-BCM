import tools as t
import numpy as np
import matplotlib.pyplot as plt
import timeit

# =============================================================================
# fileraw = "E:/more ADC data/5_21_2019_1497MA_-19dB_bypass_ADC2_attenSet=000000.bin"
# filevoltage = "E:/more ADC data/5_21_2019_1497MA_-19dB_bypass_ADC2_attenSet=000000.npy"
# # fileavg = '../ADC_DATA/8_3/1400_3000_100_avg_amplitude.npy'
# 
# 
# fs=3e9
# fo = 1.497e9
# foff = -38720 # 10 Hz accuracy (calculated with 100ms fft)
# int_time = .5*1e-3
# fsr= 1.35 # ADC32RF45 fullscale range (volts)
# ncalc = 400
# filter_ntaps = 32
# npt2n = 2**29 #= 536870912, this value is around 180ms around 3GHz sampling.
# title = ""
# ms = 1480e6 #mixer frequency for DDC calcs
# 
# #DDC settings.
# DDC = False
# bits=12
# band = "single"
# decim=12
# if DDC:
#     bits=16
#     dtype = np.uint16
#     fs=fs/decim
#     nchannels=4
# else:
#     dtype = np.int16
#     nchannels=2
# 
# #open files
# t.read_binary(infile=fileraw, outfile=filevoltage, bits=bits, fsr=1.35, send='mid', dtype=dtype, nchannels=nchannels)
# 
# if DDC==False:
#     A000000, B000000 = t.open_binary(filevoltage, nchannels=nchannels)
# else:
#     AI, AQ, BI, BQ = t.open_binary(filevoltage, nchannels=nchannels)
# 
# =============================================================================

fs = 3e9
fo = 1.425e9
int_time = 5e-4
run_time = 300*int_time
amp = .4
phi = 0
offset = np.pi/16
bits = 12
#noise = [-70,-80,-90,-100,-110]
noise = -110
IQ = True
raw = True
fsr = 1.35

# =============================================================================
# phaseAvgs = []
# phaseSigmas = []
# binnedAvgs = []
# binnedSigmas = []
# hypotAvgs = []
# hypotSigmas = []
# 
# for n in noise:
#         
#     print('Running analysis for {:} dB white noise'.format(n))
#     
#     cat = " at {:} dB of white noise, {:} bit signal resolution, {:} volt amplitude".format(n, bits, amp)
#     
#     #generate data
#     AI, AQ, BI, BQ = t.signal_gen(fs, fo, amp, phi, offset, bits, n, run_time, int_time, IQ)
#     
#     #compare amplitude reconstruction
#     phaseAvgA, phaseAvgB = t.phase_rotation(fsr*AI/2**(bits-1), fsr*AQ/2**(bits-1), 
#                                             fsr*BI/2**(bits-1), fsr*BQ/2**(bits-1), 
#                                             fs, fo, int_time, not raw)
#     
#     truth = amp*np.ones(len(phaseAvgA))
#     phaseSigma, phaseAsym = t.amplitude_correlation(phaseAvgA, phaseAvgB,
#                                                        title='Phase rotation self correlation'+cat,
#                                                        xlabel='Channel A',
#                                                        ylabel='Truth',
#                                                        hist=True, scatter = False)
#     phaseAvgs.append([phaseAvgA, phaseAvgB])
#     phaseSigmas.append(phaseSigma)
#     
#     binnedAvgA, binnedAvgB = t.binned_phase_rotation(fsr*AI/2**(bits-1), fsr*AQ/2**(bits-1),
#                                                      fsr*BI/2**(bits-1), fsr*BQ/2**(bits-1),
#                                                      fs, int_time, not raw)
#     binnedSigma, binnedAsym = t.amplitude_correlation(binnedAvgA, binnedAvgB,
#                                                           title='Pulse-wise phase rotation self correlation'+cat,
#                                                        xlabel='Channel A',
#                                                        ylabel='Truth',
#                                                        hist=True, scatter = False)
#     binnedAvgs.append([binnedAvgA, binnedAvgB])
#     binnedSigmas.append(binnedSigma)
#     
#     hypotAvgA, hypotAvgB = t.hypot_recon(fsr*AI/2**(bits-1), fsr*AQ/2**(bits-1),
#                                             fsr*BI/2**(bits-1), fsr*BQ/2**(bits-1),
#                                             fs, int_time, not raw)
#     hypotSigma, hypotAsym = t.amplitude_correlation(hypotAvgA, hypotAvgB,
#                                                        title='Hypotenuse reconstruction self correlation'+cat,
#                                                        xlabel='Channel A',
#                                                        ylabel='Truth',
#                                                        hist=True, scatter = False)
#     hypotAvgs.append([hypotAvgA, hypotAvgB])
#     hypotSigmas.append(hypotSigma)
# =============================================================================

phaseTimes=[]
hypotTimes=[]
lengths=np.geomspace(1e2,1e7,50)
loops=10

count=1
for l in lengths :
    print("loop {}".format(count))
    count+=1
    
    setup = """import numpy as np; import tools as t; l={}; fs = 3e9; fo = 1.425e9; fsr=1.35;
int_time = 5e-4; amp = .4; phi = 0;offset = np.pi/16; bits = 12; noise = -110; IQ = True; raw = True;
AI, AQ, BI, BQ = t.signal_gen(fs, fo, amp, phi, offset, bits, noise, l/fs, int_time, IQ)
    """.format(l)

    #AI, AQ, BI, BQ = t.signal_gen(fs, fo, amp, phi, offset, bits, noise, l/fs, int_time, IQ)
    
    
    p = """phaseAvgA, phaseAvgB = t.phase_rotation(fsr*AI/2**(bits-1), fsr*AQ/2**(bits-1),
    fsr*BI/2**(bits-1), fsr*BQ/2**(bits-1), fs, fo, int_time, not raw)"""
    
    h = """hypotAvgA, hypotAvgB = t.hypot_recon(fsr*AI/2**(bits-1), fsr*AQ/2**(bits-1),
    fsr*BI/2**(bits-1), fsr*BQ/2**(bits-1), fs, int_time, not raw)"""
    
    pTime = timeit.timeit(stmt=p, setup=setup, number=loops)/(2*loops)
    phaseTimes.append(pTime)
    hTime = timeit.timeit(stmt=h, setup=setup, number=loops)/(2*loops)
    hypotTimes.append(hTime)
    
a,b = t.power_fitter(lengths, np.array(phaseTimes), title="Power Law Growth fit", dataLabel="Phase reconstruction")
a,b = t.power_fitter(lengths, np.array(hypotTimes), title="Power Law Growth fit",dataLabel="V_rms reconstruction")
a,b = t.exp_fitter(lengths, np.array(phaseTimes), title="Exponential growth fit",dataLabel="Phase reconstruction")
a,b = t.exp_fitter(lengths, np.array(hypotTimes), title="Exponential growth fit", dataLabel="V_rms reconstruction")
a,b = t.quasiPoly_fitter(lengths, np.array(phaseTimes), title="Quasi-polynomial growth fit",dataLabel="Phase reconstruction")
a,b = t.quasiPoly_fitter(lengths, np.array(hypotTimes), title="Quasi-polynomial growth fit", dataLabel="V_rms reconstruction")