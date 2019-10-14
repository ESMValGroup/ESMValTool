# #############################################################################
# rainfarm.jl
# Authors:       J. von Hardenberg (ISAC-CNR, Italy)
#                E. Arnone (ISAC-CNR, Italy)
# #############################################################################
# Description
# ESMValTool diagnostic calling the RainFARM library written in Julia
# (by von Hardenberg, ISAC-CNR, Italy).
# RainFARM is a stochastic precipitation downscaling method, further adapted
# for climate downscaling.
#
# Required
# CDO
# Julia language: https://julialang.org
# RainFARM Julia library: https://github.com/jhardenberg/RainFARM.jl
#
# Optional
#
# Caveats
#
# Modification history
#    20190810-vonhardenberg_jost: rewritten in pure Julia, no R
#    20181210-vonhardenberg_jost: cleanup and using juliacall
#    20180508-arnone_enrico: Conversion to v2.0
#    20170908-arnone_enrico: 1st github version
#
# ############################################################################

import YAML
using RainFARM
using Printf

function provenance_record(infile)
    xprov = Dict( "ancestors" => infile,
        "authors" => ["vonhardenberg_jost", "arnone_enrico"],
        "references" => ["donofrio14jh", "rebora06jhm",
                         "terzago18nhess"],
        "projects" => ["c3s-magic"],
        "caption" => "RainFARM precipitation downscaling",
        "statistics" => ["other"],
        "realms" => ["atmos"],
        "themes" => ["phys"],
        "domains" => ["reg"]
    )
    return(xprov)
end

let
diag_scripts_dir = ENV["diag_scripts"]
include(joinpath(diag_scripts_dir, "shared/external.jl"))

settings = YAML.load_file(ARGS[1])

metadata = YAML.load_file(settings["input_files"][1])
climofiles = collect(keys(metadata))
climolist = metadata[climofiles[1]]
varname = climolist["short_name"]
diag_base = climolist["diagnostic"]

println(diag_base, ": starting routine")
println(diag_base, ": creating work and plot directories")
work_dir = settings["work_dir"]
run_dir = settings["run_dir"]
mkpath(work_dir)
mkpath(run_dir)
cd(run_dir)

# setup provenance file and list
provenance_file = joinpath(run_dir, "diagnostic_provenance.yml")
provenance = Dict()

# Reading parameters from the settings
nf = get(settings, "nf", 2)
slope = get(settings, "slope", 0.0)
kmin = get(settings, "kmin", 1)
nens = get(settings, "nens", 1)
weights_climo = get(settings, "weights_climo", "")
conserv_glob = get(settings, "conserv_glob", false)
conserv_smooth = get(settings, "conserv_smooth", true) 
auxiliary_data_dir = get(settings, "auxiliary_data_dir", "") 

if weights_climo isa Bool # Compatibility with old standard
    weights_climo = ""
end

# Conservation options
if (conserv_glob)
    println("Conserving global field")
elseif (conserv_smooth)
    println("Smooth conservation")
else
    println("Box conservation")
end

for (infile, value) in metadata
    (infilename, ) = splitext(basename(infile))
    outfilename = joinpath(work_dir, infilename * "_downscaled")

    println(diag_base, ": calling RainFARM for ", infilename)

    (pr, lon_mat, lat_mat) = read_netcdf2d(infile, varname)

    # Ensure grid is square and with even dims
    nmin = min(size(pr)[1], size(pr)[2])
    nmin = floor(Int, nmin / 2) * 2
    pr = pr[1:nmin, 1:nmin, :]
    if (ndims(lon_mat) == 1)
        lon_mat = lon_mat[1:nmin]
        lat_mat = lat_mat[1:nmin]
    else
        lon_mat = lon_mat[1:nmin, 1:nmin]
        lat_mat = lat_mat[1:nmin, 1:nmin]
    end

    (lon_f, lat_f) = lon_lat_fine(lon_mat, lat_mat, nf);

    # Automatic spectral slope
    if (slope == 0.)
        (fxp, ftp)=fft3d(pr)
        slope =fitslopex(fxp, kmin=kmin)
        println("Computed spatial spectral slope: ", slope)
    else
        println("Fixed spatial spectral slope: ", slope)
    end

    if weights_climo != ""
        if weights_climo[1] != '/'
            weights_climo = joinpath(auxiliary_data_dir, weights_climo)
        end
        println("Using external climatology for weights: ", weights_climo)
        fileweights = joinpath(work_dir, infilename * "_w.nc")

        ww = rfweights(weights_climo, infile, nf,
                       weightsfn = fileweights, varname = varname,
                       fsmooth = conserv_smooth)
    else
       ww = 1.
    end

    for iens=1:nens
        println("Realization ", iens)
        rd=rainfarm(pr, slope, nf, ww, fglob = conserv_glob,
                    fsmooth = conserv_smooth, verbose=true)
        fname=@sprintf("%s_%04d.nc", outfilename, iens);
        write_netcdf2d(fname, rd, lon_f, lat_f, varname, infile)
        xprov = provenance_record(infile)
        provenance[fname] = xprov
     end
end

# Write provenance file
create_yaml(provenance, provenance_file)
end
