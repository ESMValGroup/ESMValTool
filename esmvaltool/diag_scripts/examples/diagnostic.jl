# #############################################################################
# diagnostic.jl
# Authors: J. von Hardenberg (ISAC-CNR, Italy)
# #############################################################################
# Description
# Example of ESMValTool diagnostic written in Julia
#
# Modification history
#    20190807-hard_jo written for v2.0
# ############################################################################

import YAML
import JSON

function provenance_record(infile)
    xprov = Dict("ancestors" => infile,
                 "authors" => ["hard_jo", "arno_en"],
                 "references" => ["zhang-2011"],
                 "projects" => ["crescendo", "c3s-magic"],
                 "caption" => "Example diagnostic in Julia",
                 "statistics" => ["other"],
                 "realms" => ["atmos"],
                 "themes" => ["phys"],
                 "domains" => ["global"])
    return(xprov)
end

function compute_diagnostic(metadata, diag_base, parameter, work_dir)
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

        # Fake output file
        outfile = string(work_dir, "/", "outfile.txt")
        open(outfile, "w") do io
            println(io, "diagnostic output")
        end 
        # Create provenance record for the output file
        xprov = provenance_record(infile)
        provenance[outfile] = xprov
    end
    return provenance
end

# This is useful to avoid the "global scope debacle"
# https://github.com/JuliaLang/julia/issues/28789
# i.e. the following scope is local, as if we were in a function, not a script
let

diag_scripts_dir = ENV["diag_scripts"]
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

# Reading an example parameter from the settings
parameter = settings["parameter1"]

# Compute the main diagnostic
provenance = compute_diagnostic(metadata, diag_base, parameter, work_dir)

# setup provenance file
provenance_file = string(run_dir, "/diagnostic_provenance.yml")

# Write provenance file
open(provenance_file, "w") do io
    JSON.print(io, provenance, 4)
end

end
