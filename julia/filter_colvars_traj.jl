#!/usr/bin/env julia
using Printf

open("gawtm_eabf.colvars.traj") do fInput
  open("gawtm_eabf_production.traj", "w") do fOutput
    for line in eachline(fInput)
      line = strip(line)
      fields = split(line, keepempty=false)
      if length(fields) > 0 && fields[1] != "#"
        step = parse(UInt64, fields[1])
        if step > 6000000
          println(fOutput, line)
        end
      end
    end
  end
end
