
module SubModules

using DSP, Interpolations, Dierckx, Statistics, LinearAlgebra;

# ============================================================================ #
  export freqcnt
  """
      This function calculates the frequency content of input traces.

    (z, f) = freqcnt(x, dt);

      Inputs:
        ∙ x ....... input traces (columns are traces)
        ∙ dt ...... sampling interval (in s)

      Outputs:
        ∙ z ....... complex frequency response (fft) of input signals
                    (columns are traces)
        ∙ f ....... frequency vector (in Hz)

  """
  function freqcnt(x, dt)
    (n, m) = (size(x, 1), size(x, 2));
    if ( (n == 1) | (m == 1) )
      x = x[:];
      (n, m) = (size(x, 1), size(x, 2));
    end;
    N  = nextpow(2, n);

    # zero padding:
    if (N != m)
      x = [x; zeros(N-n, m)];
    end;

    df   = 1 / (N * dt);    # frequency sampling interval (Hz).
    fnyq = 1 / (2dt);       # Nyquist frequency (Hz).
    f    = collect(-fnyq: df: fnyq-df)[:];  # frequency vector (Hz).

    z = complex(zeros(N, m));
    z = DSP.fftshift(DSP.fft(x, 1), 1);

    i0 = round(Int, N/2) + 1;
    f = f[i0: end][:];
    z = z[i0: end, 1: m];

    # return output:
    return z, f;
  end

# ============================================================================ #
  export uniqueidx
  """
      This function was originally written by: Bogumił Kamiński, which returns
      the indices of unique elemetns of an input array.

    idxs = uniqueidx(x);

      Input:
        ∙ x ....... input 1D array.
      Output:
        ∙ idxs .... indices of unique elements.

      It is equal to, but faster than, the following function:
      ∙ findfirst.(isequal.(unique(x)), [x])

  """
  function uniqueidx(x::AbstractArray{T}) where T
    x = x[:];
    uniqueset = Set{T}();
    ex = eachindex(x);
    idxs = Vector{eltype(ex)}();
    for i in ex
      xi = x[i];
      if !(xi in uniqueset)
        push!(idxs, i);
        push!(uniqueset, xi);
      end;
    end;

    return idxs;
  end

# ============================================================================ #
  export snrFcn
  """
      This function calculates a continuous function of signal-to-noise rations
      using the pseudo energy of input signal.

    snr = snrFcn(din, dt, wb, wa);

      Inputs:
        din ......... input data (columns are traces)
        dt .......... sampling interval (s)
        wb .......... length of the preceeding window (s)
        wa .......... length of the post window (s)

      Output:
        snr ......... snr functions (columns are snr traces)

  """
  function snrFcn(din, dt, wb, wa)
    (n, m)   = (size(din, 1), size(din, 2));
    snr      = zeros(n, m);
    (nb, na) = (round(Int, wb/dt), round(Int, wa/dt));
    N = n - (na+nb) + 2;

    (ib, ia) = (zeros(nb, N), zeros(na, N));
    (ib, ia) = (repeat((1: nb), 1, N), repeat((1: na), 1, N));
    ib = ib + repeat(reshape((0: N-1), 1, N) , nb, 1);
    ia = ia + repeat(reshape((0: N-1), 1, N) , na, 1) .+ (nb - 1);

    (exb, exa) = (din[ib, 1: m].^2, din[ia, 1: m].^2);
    snr[nb: end-na+1, 1: m] =
      20.0 .* log10.(mean(exa, dims = 1) ./ mean(exb, dims = 1));

    # output:
    return snr;
  end

# ============================================================================ #
  export fftshift
  """
      This function apply time shift(s) to input trace(s) using fourier transform.

    dout = fftshift(din, dt, δt);

      Inputs:
        din ....... input data (columns are traces)
        dt ........ sampling interval (s)
        δt ........ time advance/delay (s)
    # NOTE: if δt is a scalar then same time shift will be applied to all traces.
    # NOTE: Otherwise, length(δt) has to be equal to the number of input traces.

      Output:
        dout ...... time-shifted signal.

  """
  function fftshift(din, dt, δt)
    npt = size(din, 1);
    npow = nextpow(2, npt);
    if (npow > npt)
      din = [din; zeros(npow-npt, size(din, 2))];
    end;
    n = size(din, 1);

    df = 1 / (n * dt);
    fnyq = 1 / 2dt;
    f = [collect(0.: df: fnyq-df)[:]; collect(-fnyq: df: -df)[:]];

    # fourier trnasform:
    z = DSP.fft(din, 1);

    # apply time shift:
    if (length(δt) == 1)
      fv = f[:] * ones(1, size(din, 2));
      zs = exp.(-2pi * im * fv * δt) .* z;

    elseif (length(δt) == size(din, 2))
      δtv = ones(n, 1) * reshape(δt, 1, length(δt));
      fv  = f[:] * ones(1, size(din, 2));
      zs = exp.(-2pi * im * (fv .* δtv)) .* z;
    elseif length(δt) != size(din, 2)
      error(
        "δt has to be a scalar or length(δt) must be equal to the number of traces!"
      );
    end;

    # inverse fourier transform:
    dout = real.(DSP.ifft(zs, 1))[1: npt, 1: size(din, 2)];

    # return time-shifted traces:
    return dout;
  end

# ============================================================================ #
  export bwbp
  """
      This function uses DSP.jl to apply the Butterworth bandpass filter.

      dfilt = bwbp(din, fs, flow, fhigh, npole, r);

    Inputs:\n
      din ......... input traces (columns are traces)
      fs .......... sampling frequency (Hz)
      flow ........ minimum corner frequency (Hz)
      fhigh ....... maximum corner frequency (Hz)
      npole ....... number of poles
      r ........... ratio of cosine taper (between 0 and 1)

    Output:\n
      dfilt ....... filtered traces (columns are traces)

  """
  function bwbp(din, fs, flow, fhigh, npole, r)

    if ( (r < 0.0) | (r > 1.0) )
      @warn(
        "
        'r' has to be between 0 and 1!
        setting the ratio of tapered section to r = 0.
        "
      );
      return
    end;
    if (flow < 0.)
      @warn(
        "
        setting the minimum corner frequency to flow = 1e-3 Hz!
        "
      );
    end;
    if (fhigh > fs/2)
      @warn(
        "
        setting the maximum corner frequency to fhigh = fs/2 - 1e-3 Hz!
        "
      );
    end;
    (n, m) = (size(din, 1), size(din, 2));

    tpr   = DSP.tukey(n, r);  # cosine taper.
    dfilt = zeros(n, m);
    dfilt = deepcopy(din) - (ones(n, 1) * mean(din, dims = 1));  # remove mean.
    dfilt = dfilt .* (tpr[:] * ones(1, m));  # apply cosine taper.
    # bandpass response:
    bpres = DSP.Bandpass(
              max(1e-3, flow), min(fhigh, fs/2-1e-3), fs = 1fs
            );
    # design the n-pole Butterworth filter:
    bwdes = DSP.Butterworth(npole);
    # filtered traces:
    dfilt = DSP.Filters.filtfilt(
              DSP.Filters.digitalfilter(bpres, bwdes),
              dfilt
            );

    return dfilt;
  end

# ============================================================================ #
  export resample
  """
    This script uses Interpolations.jl to re-sample input traces.

      dout = resample(din, ifs, ofs);

    Inputs:\n
      din .......... input traces (columns are traces)
      ifs .......... original sampling frequency (Hz)
      ofs .......... sampling frequency of output data (Hz)

    Output:
      dout ........ resampled traces (columns are resampled traces)

  """
  function resample(din, ifs, ofs)

    (n, m) = (size(din, 1), size(din, 2));
    # input and output time axes (in sec):
    itax = (0: n-1)[:] / ifs;
    otax = collect(0: 1/ofs: itax[end])[:];

    dout = zeros(length(otax), m);
    for i = 1: m
      IntFcn  = Interpolations.interpolate(
        (itax,), din[:, i],
        Interpolations.Gridded(Interpolations.Linear())
      );
      dout[:, i] = deepcopy(IntFcn(otax)[:]);
    end;

    return dout
  end

# ============================================================================ #
  export smooth1d
  """
      This script is written to smooth the input signals using a moving window.

      dout = smooth1d(din, N, wtype);

    Inputs:\n
      din ........ input signals (columns are traces).
      N .......... length (number of samples) of the moving average window.
      wtype ...... type of the averaging window:
      "rect"    -> rectangular window,
      "hamming" -> hamming window,
      "cosine"  -> cosine window,
      "triang"  -> triangular window,
      "exp"     -> exponential window.

    Output:\n
      dout ....... output smoothed signal.

  """
  function smooth1d(din, N, wtype)

    N = round.(Int, N);
    if (mod(N, 2) == 0.0)
      N = N + 1;
    end;
    N2 = floor.(Int, N/2);

    window = zeros(N, 1);
    if (wtype == "rect")
      window = DSP.rect(N);
    elseif (wtype == "hamming")
      window = DSP.hamming(N);
    elseif (wtype == "cosine")
      window = DSP.cosine(N);
    elseif (wtype == "triang")
      window = DSP.triang(N);
    elseif (wtype == "exp")
      window = collect(0: N-1) ./ collect(1: N);
    else
      error("Please select a proper window type from the list:
      'rect', 'hamming', 'cosine', 'triang', or 'exp'.   \n");
    end;

    dout = zeros(size(din));
    for ic = 1: size(din, 2)
      dconv = conv(din[:, ic], window[:]) / sum(window);
      dout[:, ic] = dconv[1+N2: end-N2][:];
    end;

    return dout
  end

# ============================================================================ #
  export smooth2d
  """
      This script is written to smooth the input signals using a moving window.

      dout = smooth1d(din, nrow, ncol, wtype, upsample);

    Inputs:\n
      din ........... input matrix (2D arrray).
      nrow .......... number of rows used for smoothing.
      ncol .......... number of columns used for smoothing.
      wtype ...... type of the averaging window:
        "rect"    -> rectangular window,
        "hamming" -> hamming window,
        "cosine"  -> cosine window,
        "triang"  -> triangular window,
        "exp"     -> exponential window.
      upsample ....... if true, it will perform upsampling first;
                       else it has to be the upsampling factor (an integer, e.g., 3).

    Output:\n
      dout .......... smoothed output.

  """
  function smooth2d(din, nrow, ncol, wtype, upsample)

    nrow = round(Int, nrow);
    ncol = round(Int, ncol);
    if (mod(nrow, 2) == 0.)
      nrow = nrow + 1;
    end;
    if (mod(ncol, 2) == 0.)
      ncol = ncol + 1;
    end;
    hnrow = floor(Int, nrow/2);
    hncol = floor(Int, ncol/2);

    winrow = zeros(nrow, 1);
    wincol = zeros(ncol, 1);
    if (wtype == "rect")
      winrow = DSP.rect(nrow);
      wincol = DSP.rect(ncol);
    elseif (wtype == "hamming")
      winrow = DSP.hamming(nrow);
      wincol = DSP.hamming(ncol);
    elseif (wtype == "cosine")
      winrow = DSP.cosine(nrow);
      wincol = DSP.cosine(ncol);
    elseif (wtype == "triang")
      winrow = DSP.triang(nrow);
      wincol = DSP.triang(ncol);
    elseif (wtype == "exp")
      winrow = collect(0: nrow-1) ./ collect(1: nrow);
      wincol = collect(0: ncol-1) ./ collect(1: ncol);
    else
      error("Please select a proper window type from the list:
      'rect', 'hamming', 'cosine', 'triang', or 'exp'.   \n");
    end;

    if (upsample == false)
      dout  = zeros(size(din));
      dconv = conv2(winrow, wincol, din) / (sum(winrow) * sum(wincol));
      dout  = dconv[hnrow+1: end-hnrow, hncol+1: end-hncol];
      dout  = dout / maximum(abs, dout) * maximum(abs, din);

    else
      N = round(Int, upsample);
      rows = collect(1: size(din, 1));
      cols = collect(1: size(din, 2));
      spl  = Dierckx.Spline2D(rows, cols, din);

      row_fine = collect(linspace(rows[1], rows[end], N*length(rows)));
      col_fine = collect(linspace(cols[1], cols[end], N*length(cols)));
      rows_ = reshape(row_fine, length(row_fine), 1) * ones(1, length(col_fine));
      cols_ = ones(length(row_fine), 1) * reshape(col_fine, 1, length(col_fine));
      dout  = Dierckx.evaluate(spl, rows_[:], cols_[:]);
      dout  = reshape(dout[:], size(rows_, 1), size(cols_, 2));

      dconv = conv2(winrow, wincol, dout) / (sum(winrow) * sum(wincol));
      dout  = dconv[hnrow+1: end-hnrow, hncol+1: end-hncol];
      dout  = dout / maximum(abs, dout) * maximum(abs, din);

    end;

    return dout
  end

# ============================================================================ #
  export localminima
  """
      This function will return the indices to local minima of the
      input time-seris.

      imin = localminima(din, N);

    Inputs:\n
      din ...... input 1D signal.
      N ........ number of samples used to calculate local minima.

    Output:\n
      imax ..... indices of local minima.

  """
  function localminima(din, N)

    if (mod(N, 2) == 0)
      N = N + 1;
    end;
    N2 = floor(Int, N/2);
    if (N2 < 2)
      N2 = 2;
    end;

    i   = collect(1: length(din))[:];
    ∂din = [0.; diff(din)[:]];
    # sign of the slope:
    sgn  = ones(size(∂din));
    sgn[findall(∂din .<= 0.)] .= -1.;
    i = findall(diff(sgn) .== +2.);
    i = i[findall( (i .> N2) .& (i .< length(din)-N2) )];
    imat = zeros(length(i), 2N2+1);
    imat = reshape(i.-N2, length(i), 1) * ones(1, 2N2+1);
    imat = imat .+ (ones(length(i), 1) * reshape(collect(0: 2N2), 1, 2N2+1));
    imat = convert.(Int64, imat);

    Δmat  = (din[i] * ones(1, 2N2+1)) .- din[imat];
    Δbool = trues(size(Δmat));
    Δbool[findall(sign.(Δmat) .> 0.)] .= false;
    Δbool = prod(Δbool, dims = 2);
    imin  = i[findall(Δbool .== true)];

    return imin;
  end

# ============================================================================ #
  export localmaxima
  """
      This function will return the indices to local maxima of the
      input time-seris.

      imax = localmaxima(din, N);

    Inputs:\n
      din ...... input 1D signal.
      N ........ number of samples used to calculate local maxima.

    Output:\n
      imax ..... indices of local maxima.

  """
  function localmaxima(din, N)

    if (mod(N, 2) == 0)
      N = N + 1;
    end;
    N2 = floor(Int, N/2);
    if (N2 < 2)
      N2 = 2;
    end;

    i   = collect(1: length(din))[:];
    ∂din = [0.; diff(din)[:]];
    # sign of the slope:
    sgn  = ones(size(∂din));
    sgn[findall(∂din .<= 0.)] .= -1.;
    i = findall(diff(sgn) .== -2.);
    i = i[findall( (i .> N2) .& (i .< length(din)-N2) )];
    imat = zeros(length(i), 2N2+1);
    imat = reshape(i.-N2, length(i), 1) * ones(1, 2N2+1);
    imat = imat .+ (ones(length(i), 1) * reshape(collect(0: 2N2), 1, 2N2+1));
    imat = convert.(Int64, imat);

    Δmat  = (din[i] * ones(1, 2N2+1)) .- din[imat];
    Δbool = trues(size(Δmat));
    Δbool[findall(sign.(Δmat) .< 0.)] .= false;
    Δbool = prod(Δbool, dims = 2);
    imax  = i[findall(Δbool .== true)];

    return imax
  end

# ============================================================================ #
  export stalta
  """
      This function calculates the modified STA/LTA ratio.

      sl = stalta(data, dt, ws, wl, method);

    Inputs: \n
      data ..... input data (columns are traces).
      dt ....... sampling interval (in s).
      ws ....... length of the short-term averaging (STA) woindow (in s).
      wl ....... length of the long-term averaging (LTA) window (in s).
      method ... method used to calculate STALTA: "Allen78" or "envelope".

    Output: \n
      sl ....... STA/LTA functions.

    * Ramin M. H. Dokht (2018)

  """
  function stalta(data, dt, ws, wl, method)

  if (wl < ws)
    error("wl has to be less than ws!");
  end;

  # size of the input data:
  (n, m) = (size(data, 1), size(data, 2));

  ns = round.(Int, ws/dt);
  nl = round.(Int, wl/dt);

  sl    = zeros(n, m);
  for i = nl+1: (n + 1 - ns)
    # short-term window:
    isb = i;
    ise = i + ns - 1;
    sta = data[isb: ise, 1: m];

    # long-term window:
    ilb = i - nl + 1;
    ile = i;
    lta = data[ilb: ile, 1: m];

    if (method == "Allen78")
      dsta = data[isb-1: ise, 1: m];
      dlta = data[ilb-1: ile, 1: m];

      # NOTE: weighting factor that balances two terms of final CF:
      # one term is related to the signal energy and the second term is
      # related to the signal frequency:
      ksta = cumsum(abs.(sta), dims = 1) ./ cumsum(abs.(diff(dsta, dims = 1)), dims = 1);
      klta = cumsum(abs.(lta), dims = 1) ./ cumsum(abs.(diff(dlta, dims = 1)), dims = 1);

      # characteristic functions:
      cf_sta = sta.^2 + ksta .* (diff(dsta, dims = 1)).^2;
      cf_lta = lta.^2 + klta .* (diff(dlta, dims = 1)).^2;

    elseif (method == "envelope")
      cf_sta = abs.(DSP.hilbert(sta));
      cf_lta = abs.(DSP.hilbert(lta));

    end;

    sl[i: i, 1: m] = (sum(cf_sta, dims = 1) ./ sum(cf_lta, dims = 1)) * (nl / ns);
  end;

  return sl;
  end

# ============================================================================ #
  export kurtosisFcn
  """
      This script is written to calculate the kurtosis function of input data.

      cf = kurtosisFcn(data, dt, wa, wintype);

    Inputs: \n
      data ............ input data (columns are traces).
      dt .............. sampling interval (in s).
      wa .............. length of the averaging window (in s).
    # NOTE: it is strognly recommended that 'wa' does not exceed the P-S delay.
      wintype ..... type of the averaging window:
        * "rect" -> rectangular window,
        * "exp"  -> exponential window.

    Output: \n
      cf .......... kurtosis function using the central moment of order 4.

    see: Baillard et al. (2014), An Automatic Kurtosis-Based P- and S-Phase \n
    Picker Designed for Local Seismic Network.

    * Ramin M. H. Dokht (2018)
  """
  function kurtosisFcn(data, dt, wa, wintype)


  # size of the input data:
  (n, m) = (size(data, 1), size(data, 2));

  cf = zeros(n, m);
  na = round.(Int, wa/dt);

  # exponential or rectangular window of size na:
  if (wintype == "exp")
    win = collect(0: na-1) ./ collect(1: na);
  elseif (wintype == "rect")
    win = ones(na, 1);
  else
    error(
      "Please select the window type from the following options:
        'rect': rectangular window
        'exp': exponential window
      "
    );
  end;

  # calculate kurtosis function:
  cf = zeros(n, m);
  i  = zeros(na, n-na+1);
  i  = (collect(1: na)[:]) * ones(1, size(i, 2));
  i  = i + ( ones(size(i, 1), 1) * reshape(collect(0: size(i, 2)-1), 1, size(i, 2)) );
  i  = round.(Int, i);
  d  = deepcopy(data[i, 1: m]);
  d  = d .* ( (win[:] * ones(1, size(i, 2))) .*  ones(size(d)));
  cf[na: end, 1: m] = (sum(d.^4 , dims = 1) / (na-1)) ./ (sum(d.^2, dims = 1) / (na-1)).^2;
  cf[1: na-1, 1: m] = ones(na-1, 1) * reshape(cf[na, 1: m], 1, m);

  # return the kurtosis function:
  return cf;
  end

# ============================================================================ #
  export scale_kurtosisFcn
  """
      This script is written to calculate the kurtosis function of input data.

      scf = scale_kurtosisFcn(cf, dt, ws);

    Inputs: \n
      cf .............. input kurtosis functions (columns are traces).
      dt .............. sampling interval (in s).
      ws .............. smoothing window length (in s).

    Output: \n
      scf ......... modified/scaled kurtosis functions (amplitudes are scaled by
      pushing down the values by the amplitude of following local maxima).

    see: Baillard et al. (2014), An Automatic Kurtosis-Based P- and S-Phase \n
    Picker Designed for Local Seismic Network.

    * Ramin M. H. Dokht (2018)
  """
  function scale_kurtosisFcn(cf, dt, ws)

  # size of the input data:
  (n, m) = (size(cf, 1), size(cf, 2));
  ns = round(Int, ws / dt);

  # -----------------------------
  # clean the initial cf of all
  # strictly negative gradients:
  # -----------------------------
  ∂cf  = diff(cf, dims = 1);
  δ    = zeros(size(∂cf));
  δ[findall(∂cf .>= 0.)] .= 1.0;
  δ[findall(∂cf .< 0.)]  .= 0.0;
  cum_ = cumsum(δ.*∂cf, dims = 1);

  F2 = zeros(n, m);
  F2[1, 1: end] = deepcopy(cf[1, 1: end]);
  F2[2: end, 1: end] = (ones(n-1, 1) * reshape(F2[1, :], 1, m)) + cum_;
  F2 = SubModules.smooth1d(F2, ns, "cosine");

  # ------------------------
  # remove the linear trend:
  # ------------------------
  a  = ones(n, 1) * reshape(F2[end, 1: end] .- F2[1, 1: end], 1, m) ./ (n .- 1.0);
  b  = ones(n, 1) * reshape(F2[1, 1: end], 1, m);
  k  = reshape(collect(1: n), n, 1) * ones(1, m);
  F3 = zeros(size(F2));
  F3 = F2 .- (a .* (k .- 1.0) .+ b);
  #F3 = SubModules.smooth1d(F3, ns, "cosine");

  # ---------------------------------------
  # scaling the amplitudes by pushing down
  # all values by the amplitude of the
  # following local maxima:
  # ---------------------------------------
  scf = zeros(size(F3));
  for ic = 1: m
    # find local maxima for the smooth CF:
    imax = SubModules.localmaxima(F3[:, ic][:], ns);
    for ii = 1: length(imax)-1
      k = findall( ( (1: n) .> imax[ii] ) .& ( (1: n) .<= imax[ii+1] ) );
      scf[k, ic] = F3[k, ic] .- F3[imax[ii+1], ic];
    end;
    scf[findall(scf[:, ic] .>= 0.0), ic] .= 0.0;
  end;

  # return scaled kurtosis functions:
  return scf;
  end

# ============================================================================ #
  export pol_covmat
  """
      This function is written to calculate polarization parameters for input
      three-component data using the eigen values of the COVARIANCE MATRIX.

      (rec, dop, inc, dip, l) = pol_covmat(data, dt, wa, wintype);

    Inputs: \n
      data ......... 3-component data (columns are traces):
        * (1st column, 2nd column, 3rd column) = (east, north, vertical)

      dt ........... sampling interval (in s).
      wa ........... length of the averaging window (in s).
      wintype ...... type of the averaging window:
        * "rect" -> rectangular window,
        * "exp"  -> exponential window.

    Output: \n
      rec ........ degree of rectilinearity.
      dop ........ degree of polarization to identify body waves.
      inc ........ vertical incidence angle (in rad).
      dip ........ dip angle of the max polarization (in rad).
      l .......... the variations in the eigen values (λi) of covariance matrix.
        * (1st column, 2nd column, 3rd column) = (λ1, λ2, λ3)
        * λ1 ≤ λ2 ≤ λ3

    see: \n
    1. Baillard et al. (2014), An Automatic Kurtosis-Based P- and S-Phase  \n
    Picker Designed for Local Seismic Network. \n
    2. Kaur et al. (2013), Detection and Identification of Seismic P-waves \n
    using Artificial Neural Networks. \n

    * Ramin M. H. Dokht (2018).
  """
  function pol_covmat(data, dt, wa, wintype)

  # size of the input data:
  (n, m) = (size(data, 1), size(data, 2));

  if (m != 3)
    error("Input data have less/or more than three components.");
  end;

  # calculate number of samples in the sliding window:
  na  = round.(Int, wa/dt);
  if (mod(na, 2) == 0.0)
    na = na + 1;
  end;

  # exponential or rectangular window of size na:
  if (wintype == "exp")
    win = collect(0: na-1) ./ collect(1: na);
  elseif (wintype == "rect")
    win = ones(na, 1);
  else
    error("Please select the window type from the following options:
    'rect': rectangular window
     'exp': exponential window");
  end;

  rec = zeros(n, 1);
  dop = zeros(n, 1);
  inc = zeros(n, 1);
  dip = zeros(n, 1);
  l   = zeros(n, m);

  for i = na: n
    ib = i - na + 1;
    ie = i;
    x  = data[ib: ie, 1] .- mean(data[ib: ie, 1]);   # east component.
    x  = x .* win;
    y  = data[ib: ie, 2] .- mean(data[ib: ie, 2]);  # north component.
    y  = y .* win;
    z  = data[ib: ie, 3] .- mean(data[ib: ie, 3]);  # vertical component.
    z  = z .* win;

    # covariance matrix:
    C = zeros(3, 3);
    C = [cov(x, x) cov(x, y) cov(x, z);
         cov(y, x) cov(y, y) cov(y, z);
         cov(z, x) cov(z, y) cov(z, z)];

    # calculate the eigenvalues and eigenvectors:
    (λ, v) = LinearAlgebra.eigen(real.(C));
    # sort in descending order:
    isort  = sortperm(λ);
    λ = λ[isort[end: -1: 1]][:];
    v = v[1: end, isort[end: -1: 1]];

    # calculate the rectilinearity:
    rec[i] = 1.0 - ( (λ[2] + λ[3]) / (2.0 * λ[1]) );

    # calculate the degree of polarization:
    dop[i] = ((λ[3]-λ[2])^2 + (λ[2]-λ[1])^2 + (λ[1]-λ[3])^2) /
             (2.0 * (λ[1]+λ[2]+λ[3])^2);

    # calculate vertical incidence angle:
    inc[i] = abs.(atan(sqrt.(v[2, 1]^2 + v[1, 1]^2) / v[3, 1]));

    # dip angle of the max polarization:
    dip[i] = atan(v[3, 1] / sqrt(v[2, 1]^2 + v[1, 1]^2));

    # the largest eigen value of the covariance matrix:
    l[i, :] = [λ[1] λ[2] λ[3]];

  end;

  rec[1: na-1] .= rec[na];
  dop[1: na-1] .= dop[na];
  inc[1: na-1] .= inc[na];
  dip[1: na-1] .= dip[na];
  l[1: na-1, 1: m] = repeat(l[na: na, 1: m], na-1, 1);

  return rec, dop, inc, dip, l;
  end

# ============================================================================ #

end
