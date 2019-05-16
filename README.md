# Seismology
Tools for research in Seismology

The current repository includes the following scripts and modules:

1. get_response_SACPZs.jl
This function extracts the instrument response information from a SAC Poles and Zeros response file for a given DateTime.

2. removeIR.jl
This function deconvolve the instrument response from the input trace(s) using the SAC Poles and Zeros of instrument response.

3. SubModules.jl
This module includes the following functions:
3.1. freqcnt:
  It calculates the frequency content of input trace(s) using the DSP package.
3.2. uniqueidx:
  It returns the indices of unique elements of an input 1D array.
3.3. snrFcn:
  It calculates the SNR function(s) of input trace(s).
3.4. fftshift:
  It applies time-shift(s) to input trace(s) using the fast fourier transform and the DSP package.
3.5. bwbp:
  It applies the Butteworth bandpass filter to input trace(s) using the DSP package.
3.6. resample:
  Ressampling the input trace(s) using the Interpolations package.
3.7. smooth1d:
  Smoothing the input trace(s) using a moving window.
3.8. smooth2d:
  Smoothing 2D array using 2D convolution of input array with two moving windows.
3.9. localminima:
  It returns the indicies of local minima of an input trace.
3.10. localmaxima:
  It returns the indicies of local maxima of an input trace.
3.11. stalta:
  It calculates the STA/LTA functions of input traces using the original traces and/or their envelopes.
3.12. kurtosisFcn:
  Calculating the kurtosis functions of input traces.
3.13. scale_kurtosisFcn:
  Calculating the modified/scaled kurtosis functions from the original kurtosis functions.
3.14. pol_covmat:
  It calculates the polarization parameters (such as: rectilinearity, degree of polarization, and etc) for input three-component waveform data using tje eigen values of the covariance matrix.
  


