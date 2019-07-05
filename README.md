# Seismology
Julia Tools for research in Seismology

The current repository includes the following scripts and modules:

**1. get_response_SACPZs.jl**\
This function extracts the instrument response information from a SAC Poles and Zeros response file for a given DateTime.

**2. removeIR.jl**\
This function deconvolve the instrument response from the input trace(s) using the SAC Poles and Zeros of instrument response.

**3. SubModules.jl**\
This module includes the following functions:
- freqcnt:\
It calculates the frequency content of input trace(s) using the DSP package.
- uniqueidx:\
  It returns the indices of unique elements of an input 1D array.
- snrFcn:\
  It calculates the SNR function(s) of input trace(s).
- fftshift:\
  It applies time-shift(s) to input trace(s) using the fast fourier transform and the DSP package.
- bwbp:\
  It applies the Butteworth bandpass filter to input trace(s) using the DSP package.
- resample:\
  Ressampling the input trace(s) using the Interpolations package.
- smooth1d:\
  Smoothing the input trace(s) using a moving window.
- smooth2d:\
  Smoothing 2D array using 2D convolution of input array with two moving windows.
- localminima:\
  It returns the indicies of local minima of an input trace.
- localmaxima:\
  It returns the indicies of local maxima of an input trace.
- stalta:\
  It calculates the STA/LTA functions of input traces using the original traces and/or their envelopes.
- kurtosisFcn:\
  Calculating the kurtosis functions of input traces.
- scale_kurtosisFcn:\
  Calculating the modified/scaled kurtosis functions from the original kurtosis functions.
- pol_covmat:\
  It calculates the polarization parameters (such as: rectilinearity, degree of polarization, and etc) for input three-component waveform data using the eigen values of the covariance matrix.
  


