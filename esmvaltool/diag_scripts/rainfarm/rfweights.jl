#!/usr/bin/env julia

# rfweights create weights for RainFARM downscaling

# RainFARM 
# Stochastic downscaling following 
# D'Onofrio et al. 2014, J of Hydrometeorology 15 , 830-843 and
# Rebora et. al 2006, JHM 7, 724 
# Includes orographic corrections

# Implementation in Julia language
# Copyright (c) 2016, Jost von Hardenberg - ISAC-CNR, Italy

using RainFARM
using ArgParse
using Compat, Compat.Printf

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--nf", "-n"
            help = "Subdivisions for downscaling"
            arg_type = Int
            default = 2
        "--weights", "-w", "--weight", "--outfile", "-o", "--out"
            help = "Output weights filename"
            arg_type = AbstractString
            default = "weights.nc" 
        "--varname", "-v"
            help = "Input variable name (in orofile)"
            arg_type = AbstractString
            default = "" 
	"--pass2", "-p"
            help = "Previous rainfarm output for 2-pass weighting"
            arg_type = AbstractString
            default = "" 
	"--inweight", "-i"
            help = "Previous weights file to correct for 2-pass"
            arg_type = AbstractString
            default = "" 
        "orofile"
            help = "The input file to use for orography"
            arg_type = AbstractString
            required = true
        "reffile"
            help = "A reference file (e.g. the file to downscale)"
            arg_type = AbstractString
            required = true
        "--conv", "-c"
            action = :store_true
            help = "conserve precipitation using convolution"
    end

    s.description="Create weights for RainFARM downscaling"
    s.version="0.1"
    s.add_version=true

    return parse_args(s)
end

args = parse_commandline()
nf=args["nf"]
reffile=args["reffile"]
orofile=args["orofile"]
weightsfn=args["weights"]
varname=args["varname"]
pass1fn=args["pass2"]
inweight=args["inweight"]
fsmooth=args["conv"]

println("Creating weights from file ",orofile)

# Create a reference gridrf.nc file (same grid as rainfarm output files)
(pr,lon_mat,lat_mat)=read_netcdf2d(reffile, varname);
# Creo la griglia fine
nss=size(pr)
if (length(nss)>=3)
    pr=pr[:,:,1]
end
println(nss)
ns=nss[1];
(lon_f, lat_f)=lon_lat_fine(lon_mat, lat_mat,nf);

rr=round.(Int,rand(1)*100000)

println("Output size: ",size(lon_f))
if(varname=="")
   varname="pr"
   run(`cdo -s  setname,pr $reffile reffile_rr.nc`)
   reffile="reffile_rr.nc"
end

# The rest is done in CDO

if(inweight=="")
println("Computing weights")
write_netcdf2d("gridrf.nc",reshape(pr,ns,ns,1),lon_f,lat_f,varname,reffile)
run(`cdo -s timmean $orofile pr_orofile_$rr.nc`)
run(`cdo -s -f nc copy gridrf.nc gridrf_2_$rr.nc`)
run(`cdo -s -f nc remapbil,gridrf_2_$rr.nc pr_orofile_$rr.nc pr_remap_rr.nc`)
if(fsmooth)
  (prr,lon,lat)=read_netcdf2d("pr_remap_rr.nc","")
  ww=prr./smoothconv(prr,ns);
  write_netcdf2d(weightsfn,ww,lon_f,lat_f,varname,reffile)
  run(`rm -f pr_remap_rr.nc pr_orofile_$rr.nc gridrf.nc gridrf_2_$rr.nc reffile_rr.nc`)
else
  run(`cdo -s gridboxmean,$nf,$nf pr_remap_rr.nc pr_remap_gbm_$rr.nc`)
  run(`cdo -s remapnn,pr_remap_rr.nc pr_remap_gbm_$rr.nc pr_remap_nn_$rr.nc`)
  run(`cdo -s div pr_remap_rr.nc pr_remap_nn_$rr.nc $weightsfn`)
  run(`rm -f pr_remap_rr.nc pr_remap_gbm_$rr.nc pr_remap_nn_$rr.nc pr_orofile_$rr.nc gridrf.nc gridrf_2_$rr.nc reffile_rr.nc`)
  inweight=weightsfn
end
else
println("Correcting weights in ",inweight)
end

if(pass1fn!="")
println("2-pass weighting using: ",pass1fn)

run(`cdo -s timmean $pass1fn pr_remap_$rr.nc`)
run(`cdo -s gridboxmean,$nf,$nf pr_remap_$rr.nc pr_remap_gbm_$rr.nc`)
run(`cdo -s -f nc copy pr_remap_$rr.nc pr_remap_2_$rr.nc`)
run(`cdo -s remapnn,pr_remap_2_$rr.nc pr_remap_gbm_$rr.nc pr_remap_nn_$rr.nc`)
run(`cdo -s div pr_remap_2_$rr.nc pr_remap_nn_$rr.nc neweights_$rr.nc`)

run(`cdo -s div $inweight neweights_$rr.nc weightsn_$rr.nc`)
run(`mv weightsn_$rr.nc $weightsfn`)
run(`rm pr_remap_$rr.nc pr_remap_gbm_$rr.nc pr_remap_nn_$rr.nc neweights_$rr.nc pr_remap_2_$rr.nc`)
end

