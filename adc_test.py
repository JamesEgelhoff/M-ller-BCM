#import useful packages
import tools as t
import numpy as np
import matplotlib.pyplot as plt
import timeit

#store the addresses in memory where data will be loaded from and saved to
fileBinary = "E:/more ADC data/example.bin"
fileNumpy = "E:/more ADC data/example.npy"

#intialize useful variables
fs=3e9 #sampling frequency in Hz
fo = 1.497e9 #signal frequency in Hz
int_time = .5*1e-3 #length of pulse/integration window in seconds
fsr= 1.35 # ADC32RF45 fullscale range (volts)

#choose settings for read_binary and open_binary
#if the binary file consists of DDC data, we want to read the file with different settings
#since DDC and bypass data have different data formats

DDC = False #enter whether to read binary file is DDC data or bypassed data
decim=12 #enter the decimation factor is data was decimated
if DDC:
    bits=16 #bit precision
    dtype = np.uint16
    fs=fs/decim
    nchannels=4 #4 channels: 2 ADCs, each outputting an in-phase and quadrature signal if in DDC
else:
    bits=12 #bit precision
    dtype = np.int16
    nchannels=2 #2 channels: 2 ADCs, eash outputting a single real signal

#convert the .bin file to .npy and save
t.read_binary(infile=fileBinary, outfile=fileNumpy, bits=bits, fsr=1.35, send='mid', dtype=dtype, nchannels=nchannels)

#load the .npy file
if DDC==False:
    A, B = t.open_binary(fileNumpy, nchannels=nchannels)
else:
    AI, AQ, BI, BQ = t.open_binary(fileNumpy, nchannels=nchannels)

#Now that some data is loaded up you can perform an analysis.
#Note: read_binary is slow, so I recommend that if the binary file has already
#been converted and save to a .npy file you comment out the line with read_binary
#to speed up your analysis run time.

#example analysis:
#resolution analysis - calculates the relative amplitude asymmetry between
#successive pulse window pairs. Compares the relative difference in this value
#for channel A and channel B to determine the resolution of the signal generation
#and measurement setup in parts per million. This relative difference is stored
#for each pair of pulse windows and plotted in a histogram
sigma, mu, diffs = t.amplitude_asym_hist(A, B, fs=fs, pulse_width=int_time)