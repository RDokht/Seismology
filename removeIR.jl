using DSP, Statistics, Polynomials;
  """
      This program is written to remove instrument response using
      SAC Poles and Zeros of instrument response.

      ● dout = removeIR(din, p, z, fs, flow, fhigh, cons);

    Inputs: \n
      din ....... input trace.
      p ......... a vector of poles (complex values).
      z ......... a vector of zeros (complex values).
      fs ........ sampling frequency (in Hz).
      flow ...... min corner freq (Hz) used for bandpass filtering.
      fhigh ..... max corner freq (Hz) used for bandpass filtering.
      cons ...... constant = sensitivity * normalization factor (A0).
      R ......... ratio of tapered section to total length of signal (0 to 1)

    Output: \n
      dout ...... Instrument response-corrected seismogram.

    This script runs on Julia 1.0

    author: Ramin M.H. Dokht (2018).

    The MATLAB version of this code has been written by: Martin Mityska

  """
function removeIR(din, p, z, fs, flow, fhigh, cons, R)

if ( (R < 0) | (R > 1) )
  error("R has to be between 0.0 and 1.0!");
  return Array{Float64}(0);
end;

# remove mean:
din  = din[:];
din  = din[:] .- Statistics.mean(din);

# remove linear trend:
linear_fcn = Polynomials.polyfit((1: length(din))[:], din, 1);
din = din .- linear_fcn((1: length(din))[:]);
din = din .- Statistics.mean(din);
din = din .* DSP.tukey(length(din), R);

N  = length(din);
zb = round(Int, 0.5N);
za = deepcopy(zb);
if (mod(zb+N+za, 2) != 0)
  za += 1;
end;

dinpad = [
  zeros(zb, 1);
  din[:];
  zeros(za, 1)
];

N  = length(dinpad);

# Construct numerator and denominator of the transfer fcn:
df   = fs / N;
fnyq = fs / 2.0;
ws   = collect(0: N/2)[:] * df * 2pi;

tf = Dict(
  :z => complex.(ones(length(ws), 1)),
  :p => complex.(ones(length(ws), 1))
);

for iz = 1: length(z)
  tf[:z] = tf[:z] .* (ws * im .- z[iz]);
end;
for ip = 1: length(p)
  tf[:p] = tf[:p] .* (ws * im .- p[ip]);
end;

H = 1 ./ ( (tf[:z] ./ tf[:p])[:] .* cons );
dinspec = DSP.fft(dinpad[:]);
CorSpec = dinspec[1: round(Int, N/2+1)] .* H;

inan = findall(isnan.(CorSpec));
CorSpec[inan] .= complex(0.0);
CorSpec[1]     = complex(0.0);

DOUT = [CorSpec[:]; conj(CorSpec[end-1: -1: 2])[:]];
dout = real(DSP.ifft(DOUT));

# apply 4th order [zero-phase] Butterworth bandpass filter:
flow  = max(1e-4, flow);
fhigh = min(fs/2-1e-4, fhigh);
responsetype = DSP.Bandpass(flow, fhigh, fs = 1.0fs);
designmethod = DSP.Butterworth(4);
dout = DSP.Filters.filtfilt(
  DSP.Filters.digitalfilter(responsetype, designmethod),
  dout
);
dout = dout[zb+1: end-za];
dout = dout .- Statistics.mean(dout);

# Return output (instrument response-corrected signal):
return dout;

end
