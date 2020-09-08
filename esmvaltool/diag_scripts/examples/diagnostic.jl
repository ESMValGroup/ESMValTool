# #############################################################################
# diagnostic.jl
# Authors: J. von Hardenberg (ISAC-CNR, Italy)
# #############################################################################
# Description
# Example of ESMValTool diagnostic written in Julia
#
# Modification history
#    20190807-vonhardenberg_jost written for v2.0
#    20191117-vonhardenberg_jost added more realistic writing of file and plot
# ############################################################################

import YAML
import JSON
using NetCDF
# Used to write output NetCDF file with original attributes
using RainFARM
using Statistics

using PyPlot
# Avoid plotting to screen
pygui(false)

# Provides the plotmap() function
include(joinpath(dirname(@__DIR__), "shared/external.jl"))

function provenance_record(infile)
    xprov = Dict("ancestors" => [infile],
                 "authors" => ["vonhardenberg_jost", "arnone_enrico"],
                 "references" => ["zhang11wcc"],
                 "projects" => ["crescendo", "c3s-magic"],
                 "caption" => "Example diagnostic in Julia",
                 "statistics" => ["other"],
                 "realms" => ["atmos"],
                 "themes" => ["phys"],
                 "domains" => ["global"])
    return(xprov)
end

function compute_diagnostic(metadata, varname, diag_base, parameter,
                            work_dir, plot_dir)
    provenance = Dict()
    for (infile, value) in metadata
        dataset = value["dataset"]
        reference_dataset = value["reference_dataset"]
        start_year = value["start_year"]
        end_year = value["end_year"]
        exp = value["exp"]
        ensemble = value["ensemble"]
        println(diag_base, ": working on file ", infile)
        println(diag_base, ": calling diagnostic with following parameters")
        println(dataset, " ", reference_dataset, " ", start_year, " ",
                end_year, " ", exp, " ", ensemble, parameter)
        # Call here actual diagnostic
        println(diag_base, ": I am your Julia diagnostic")

        # Read the variable, lon and lat
        var = ncread(infile, varname)
        lon = ncread(infile, "lon")
        lat = ncread(infile, "lat")

        units = ncgetatt(infile, varname, "units")

        # Compute time average and add parameter
        varm = mean(var, dims = 3) .+ parameter

        # Output filename
        outfile = string(work_dir, "/", varname, "_", dataset, "_", exp, "_",
                         ensemble, "_", start_year, "-",
                         end_year, "_timmean.nc")

        # Use the RainFARM function write_netcdf2d to write variable to
        # output file copying original attributes from infile
        write_netcdf2d(outfile, varm, lon, lat, varname, infile)

        # Create provenance record for the output file
        xprov = provenance_record(infile)

        # Plot the field
        plotfile = string(plot_dir, "/", varname, "_", dataset, "_", exp, "_",
                         ensemble, "_", start_year, "-",
                         end_year, "_timmean.png")
        title = string("Mean ", varname, " ", dataset, " ", exp, " ", ensemble,
                       " ", start_year, "-", end_year)
        plotmap(lon, lat, var, title = title, proj = "robinson", clabel = units)
        savefig(plotfile)
        xprov["plot_file"] = plotfile
        provenance[outfile] = xprov
    end
    return provenance
end

function main(settings)

    metadata = YAML.load_file(settings["input_files"][1])
    climofiles = collect(keys(metadata))
    climolist = metadata[climofiles[1]]
    varname = climolist["short_name"]
    diag_base = climolist["diagnostic"]

    println(diag_base, ": starting routine")
    println(diag_base, ": creating work and plot directories")
    work_dir = settings["work_dir"]
    plot_dir = settings["plot_dir"]
    run_dir = settings["run_dir"]
    mkpath(work_dir)
    mkpath(run_dir)
    mkpath(plot_dir)
    cd(run_dir)

    # Reading an example parameter from the settings
    parameter = settings["parameter1"]

    # Compute the main diagnostic
    provenance = compute_diagnostic(metadata, varname, diag_base,
                                    parameter, work_dir, plot_dir)

    # setup provenance file
    provenance_file = string(run_dir, "/diagnostic_provenance.yml")

    # Write provenance file
    open(provenance_file, "w") do io
        JSON.print(io, provenance, 4)
    end
end

settings = YAML.load_file(ARGS[1])
main(settings)
