using Dates
  """
      This function is written to extract the instrument information from a
      SAC Poles and Zeros file, which is generate by rdseed, for a given date.

      ● response = get_response_SACPZs(sacpzfile, reqdate);

    Inputs: \n
      sacpzfile ........ path to the SAC Poles and Zeros file.
      reqdate .......... request date and time
                         (e.g. Dates.DateTime(2018, 01, 20, 00, 00, 00));

    Output: \n
      response ......... a dictionary array holding the following keys:
       :network -> network name,
       :station -> station name,
       :channel -> channel name,
           :lat -> instrument latitude (deg),
           :lon -> instrument longitude (deg),
           :ele -> instrument elevation (m),
            :az -> instrument azimuth (deg), measured from north,
           :dip -> instrument dip (deg), measured from the vertical line,
            :fs -> sample rate (Hz),
      :constant -> instrument constant (sensitivity-times-normalization factor),
         :Poles -> Poles of the transfer function,
         :Zeros -> Zeros of the transfer function.

    This script runs on Julia 1.0

    author: Ramin M.H. Dokht (2018)

  """
function get_response_SACPZs(sacpzfile, reqdate)

fid = open(sacpzfile, "r");
cnt = readlines(fid);
close(fid);

for il = 1: length(cnt) # remove tabs (\t) and returns (\n):
  cnt[il] = replace(cnt[il], "\n" => "");
  cnt[il] = replace(cnt[il], "\t" => " ");
end;

# find matching date (if multiple responses exist):
START = Array{Dates.DateTime}(undef, 0);
END   = Array{Dates.DateTime}(undef, 0);
(bool, il) = (false, Array{Int64}(undef, 0));
n     = 1;
while (n <= length(cnt))

  if (typeof(findfirst("START", cnt[n])) != Nothing)
    icol = findfirst(isequal(':'), cnt[n]);
    if !isempty(icol)
      icol  = icol[1] + 1;
      if (cnt[n][end-1: end] == "60")
        cnt[n] = "$(cnt[n][1: end-2])59";
      end;
      if (cnt[n+1][end-1: end] == "60")
        cnt[n+1] = "$(cnt[n+1][1: end-2])59";
      end;
      start_ = Dates.DateTime(rstrip(lstrip(cnt[n][icol: end])));
      end_   = Dates.DateTime(rstrip(lstrip(cnt[n+1][icol: end])));

      if ( (reqdate >= start_) & (reqdate <= end_) )
        START = deepcopy(start_);
        END   = deepcopy(end_);
	      (bool, il) = (true, deepcopy(n));
        break;
      end;
    end;
  end;
  n += 1;
end;

if !bool
  return Dict();

else
  network  = Array{String}(undef, 0);
  station  = Array{String}(undef, 0);
  channel  = Array{String}(undef, 0);
  lat      = Array{Float64}(undef, 0);
  lon      = Array{Float64}(undef, 0);
  ele      = Array{Float64}(undef, 0);
  az       = Array{Float64}(undef, 0);
  dip      = Array{Float64}(undef, 0);
  fs       = Array{Float64}(undef, 0);
  constant = Array{Float64}(undef, 0);
  Zeros    = Array{ComplexF64}(undef, 0);
  Poles    = Array{ComplexF64}(undef, 0);

  ib = Array{Int64}(undef, 0);
  ie = Array{Int64}(undef, 0);
  for n = il-1: -1: 1
     if (typeof(findfirst("NETWORK", cnt[n])) != Nothing)
       ib = deepcopy(n);
       break;
     end;
  end;
  for n = il+1: length(cnt)
    if (typeof(findfirst("CONSTANT", cnt[n])) != Nothing)
      ie = deepcopy(n);
      break;
    end;
  end;

  cnt = cnt[ib: ie];

  inet  = findfirst.("NETWORK", cnt);
  ista  = findfirst.("STATION", cnt);
  ich   = findfirst.("CHANNEL", cnt);
  ilat  = findfirst.("LATITUDE", cnt);
  ilon  = findfirst.("LONGITUDE", cnt);
  iele  = findfirst.("ELEVATION", cnt);
  idip  = findfirst.("DIP", cnt);
  iaz   = findfirst.("AZIMUTH", cnt);
  ifs   = findfirst.("SAMPLE RATE", cnt);
  izero = findfirst.("ZEROS", cnt);
  ipole = findfirst.("POLES", cnt);
  icons = findfirst.("CONSTANT", cnt);

  if !isempty(findall(inet .!= nothing))
    n = findall(inet .!= nothing)[1];
    network = strip(split(cnt[n], ":")[end]);
  end;
  if !isempty(findall(ista .!= nothing))
    n = findall(ista .!= nothing)[1];
    station = strip(split(cnt[n], ":")[end]);
  end;
  if !isempty(findall(ich .!= nothing))
    n = findall(ich .!= nothing)[1];
    channel = strip(split(cnt[n], ":")[end]);
  end;
  if !isempty(findall(ilat .!= nothing))
    n = findall(ilat .!= nothing)[1];
    lat = parse(Float64, strip(split(cnt[n], ":")[end]));
  end;
  if !isempty(findall(ilon .!= nothing))
    n = findall(ilon .!= nothing)[1];
    lon = parse(Float64, strip(split(cnt[n], ":")[end]));
  end;
  if !isempty(findall(iele .!= nothing))
    n = findall(iele .!= nothing)[1];
    ele = parse(Float64, strip(split(cnt[n], ":")[end]));
  end;
  if !isempty(findall(iaz .!= nothing))
    n = findall(iaz .!= nothing)[1];
    az = parse(Float64, strip(split(cnt[n], ":")[end]));
  end;
  if !isempty(findall(idip .!= nothing))
    n = findall(idip .!= nothing)[1];
    dip = parse(Float64, strip(split(cnt[n], ":")[end]));
  end;
  if !isempty(findall(ifs .!= nothing))
    n = findall(ifs .!= nothing)[1];
    fs = parse(Float64, strip(split(cnt[n], ":")[end]));
  end;
  if !isempty(findall(icons .!= nothing))
    n = findall(icons .!= nothing)[1];
    constant = parse(Float64, split(cnt[n])[end]);
  end;
  if !isempty(findall(izero .!= nothing))
    n = findall(izero .!= nothing)[1];
    nz = parse(Int64, split(cnt[n])[end]);  # no of zeros.
    for i = n+1: n+nz
      push!(Zeros,
        parse(Float64, split(cnt[i])[1]) + im*parse(Float64, split(cnt[i])[2])
      );
    end;
  end;
  if !isempty(findall(ipole .!= nothing))
    n = findall(ipole .!= nothing)[1];
    np = parse(Int64, split(cnt[n])[end]);  # no of poles.
    for i = n+1: n+np
      push!(Poles,
        parse(Float64, split(cnt[i])[1]) + im*parse(Float64, split(cnt[i])[2])
      );
    end;
  end;

  response = Dict(
    :network  => deepcopy(network),
    :station  => deepcopy(station),
    :channel  => deepcopy(channel),
    :lat      => deepcopy(lat),
    :lon      => deepcopy(lon),
    :ele      => deepcopy(ele),
    :az       => deepcopy(az),
    :dip      => deepcopy(dip),
    :fs       => deepcopy(fs),
    :constant => deepcopy(constant),
    :Zeros    => deepcopy(Zeros),
    :Poles    => deepcopy(Poles)
  );

  # outpus:
  return response;

end;

end
