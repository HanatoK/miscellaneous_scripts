#!/usr/bin/env julia
using Printf

function HistogramIndex(x, lowerboundary, upperboundary, width)
  if x < lowerboundary
    return (0, false)
  elseif x > upperboundary
    return (0, false)
  elseif x == upperboundary
    idx = Int(floor((x - lowerboundary) / width))
    return (idx, true)
  else
    idx = Int(floor((x - lowerboundary) / width)) + 1
    return (idx, true)
  end
end

lb_x = -9.5
lb_y = -9.5
width_x = 0.1
width_y = 0.1
bin_x = 190
bin_y = 190
ub_x = lb_x + width_x * bin_x
ub_y = lb_y + width_y * bin_y
histogram_x = zeros(bin_x)
histogram_xy = zeros(bin_x, bin_y)
histogram_x_meanV = zeros(bin_x)
histogram_x_meanV2 = zeros(bin_x)
pmf_x = zeros(bin_x)
beta = 1.0 / (300.0 * 0.0019872041)

open("GaMD.traj", "r") do traj
  for line in eachline(traj)
    line = strip(line)
    fields = split(line, keepempty=false)
    if length(fields) > 0 && fields[1] != "#"
      x = parse(Float64, fields[1])
      y = parse(Float64, fields[2])
      index_x, in_boundary_x = HistogramIndex(x, lb_x, ub_x, width_x)
      index_y, in_boundary_y = HistogramIndex(y, lb_y, ub_y, width_y)
      if in_boundary_x
        histogram_x[index_x] += 1
        bias_energy = parse(Float64, fields[6])
        histogram_x_meanV[index_x] += bias_energy
        histogram_x_meanV2[index_x] += bias_energy * bias_energy
        if in_boundary_y
          histogram_xy[index_x, index_y] += 1
        end
      end
    end
  end
end

open("histogram_1D.count", "w") do outfile
  for ix in 1:bin_x
    px = lb_x + (ix - 0.5) * width_x
    value = histogram_x[ix]
    @printf(outfile, "%12.7f %10d\n", px, value)
  end
end

open("histogram_1D.pmf", "w") do outfile
  for ix in 1:bin_x
    count = histogram_x[ix]
    if count > 0
      histogram_x_meanV[ix] /= count
      histogram_x_meanV2[ix] /= count
      variance = histogram_x_meanV2[ix] - histogram_x_meanV[ix] * histogram_x_meanV[ix]
      pmf_x[ix] = -1.0 * (histogram_x_meanV[ix] + 0.5 * beta * variance)
    end
  end
  max_pmf = maximum(pmf_x)
  for ix in 1:bin_x
    if histogram_x[ix] == 0
      pmf_x[ix] = max_pmf
    end
  end
  min_pmf = minimum(pmf_x)
  for ix in 1:bin_x
    px = lb_x + (ix - 0.5) * width_x
    @printf(outfile, "%12.7f %12.7f\n", px, pmf_x[ix] - min_pmf)
  end
end

open("histogram_2D.count", "w") do outfile
  for ix in 1:bin_x
    for iy in 1:bin_y
      px = lb_x + (ix - 0.5) * width_x
      py = lb_y + (iy - 0.5) * width_y
      value = histogram_xy[ix, iy]
      @printf(outfile, "%12.7f %12.7f %10d\n", px, py, value)
    end
  end
end
