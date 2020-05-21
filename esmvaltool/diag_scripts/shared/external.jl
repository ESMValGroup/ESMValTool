using PyPlot, PyCall
using NetCDF
	
"""
    create_yaml(filename, dict)
Write a dictionary to a YAML file.
"""
function create_yaml(dict::Dict{Any,Any}, filename::AbstractString)
  os = open(filename, "w")
  for (key, value) in dict
    println(os, "? ",key)
    indent = ": "
    print_yaml(os, value, indent)
  end
  close(os)
end

function print_yaml(os::IOStream, obj::Dict{String,Any}, indent::String)
  for (key, value) in obj
    print(os, indent, key, ": ")
    print_yaml(os, value, "  ")
    indent = "  "
  end
end

function print_yaml(os::IOStream, obj::String, indent::String)
  println(os, obj)
end

function print_yaml(os::IOStream, obj::Array{String}, indent::String)
  println(os)
  for i = 1:length(obj)
    println(os, indent, "- ", obj[i])
  end
end

"""
    plotmap(fname, var; ...)

Easy plotting of gridded datasets on a global map.
Plots variable `var` from a netcdf file `fname`.
Optional arguments are available to control the details of the plot.

# Arguments
- 'fname::String': netcdf filename containing the data to plot
- 'var'::String  : name of the variable to plot

Optional:
- 'lon::String': name of the lon variable ("lon")
- 'lat::String': name of the lat variable ("lat")
- 'lonb::String': name of the lon bounds variable ("lon_bnds")
- 'latb::String': name of the lat boundsvariable ("lat_bnds")
- 'title':      title string ("")
- 'cstep': 	 divisions of the colorbar axis ([])
- 'cmap':        colormap ("RdBu_r")
- 'proj':        projection. One of ["platecarree", "robinson", "mollweide"]. Defaults to "platecarree".
- 'cpad':        padding (shift) of the colorbar (0.08)
- 'sub':         matplotlib subplot option (e.g. "221" for the first panel of 2x2 subplots) ("111")
- 'clabel':      label of the colorbar (defaults to the units string read from the netcdf file)
- 'cdir':        direction of the colorbar. One of ["horizontal", "vertical"]. ("horizontal")
- 'cscale':      scaling of the colorbar (0.65)
- 'cfs':         colorbar ticks font size (12)
- 'lfs':         colorbar label font size (12)
- 'tfs':         title font size (14)
- 'tpad':        padding (shift) of the title string
- 'tweight':     weight of title font. One of [ 'normal" "bold" "heavy" "light" "ultrabold" "ultralight"]. ("normal")
- 'grid':        grid spacing (defaults to [60,30]). set to empty [] to remove gridlines.
- 'region::NTuple{4,Int64}':       region to plot in format (lon1, lon2, lat1, lat2). Defaults to global.
- 'style':       one of ["pcolormesh" "contourf"]. Defaults to "pcolormesh".
- 'levels':      contour plot levels. Can be an array or a number of levels (auto)
- 'extend':      plot levels outside `levels` range. One of ["neither", "both", "min", "max"]. Default: "neither".

Author: Jost von Hardenberg, 2019
"""
function plotmap(fname::String, var::String; lon="lon", lat="lat", 
                 lonb="lon_bnds", latb="lat_bnds", title="", cstep=[],
                 cmap="RdBu_r", proj="", cpad=0.08, tpad=24, sub=111,
                 clabel="NONE", cdir="horizontal", cscale=0.65, tfs=14,
                 cfs=12, lfs=12, tweight="normal", grid=[60,30], region=(),
                 style="pcolormesh", levels=0, extend="neither")

# pcolormesh needs cell boundaries
if style=="pcolormesh"
    try 
        lonb=ncread(fname, lonb);
        lonv=vcat(lonb[1,:],lonb[2,end])
    catch
        lonv=ncread(fname, lon);
    end
    try
        latb=ncread(fname, latb);
        latv=vcat(latb[1,:],latb[2,end])
    catch
        latv=ncread(fname, lat);
    end
else
    lonv=ncread(fname, lon);
    latv=ncread(fname, lat);
end

data=ncread(fname, var);
units=ncgetatt(fname, var, "units")
if clabel=="NONE" clabel=units end

plotmap(lonv, latv, data; title=title, cstep=cstep, cmap=cmap, proj=proj,
        cpad=cpad, tpad=tpad, sub=sub, clabel=clabel, cdir=cdir,
        cscale=cscale, tfs=tfs, cfs=cfs, lfs=lfs, tweight=tweight, grid=grid,
        region=region, style=style, levels=levels, extend=extend)

end

"""
    plotmap(lon, lat, data; ...)

Easy plotting of gridded datasets on a global map.
Plots data in 2D array `data` with longitudes `lon` and latitudes `lat`.
Optional arguments are available to control the details of the plot.

# Arguments
- 'data::Array{Float32,2}': data to plot. If a 3D array is passed, only the first frame is plotted: `data[:,:,1]`
- 'lon::Array{Float64,1}' : longitudes
- 'lat::Array{Float64,1}' : latitudes

Optional:
- 'title':      title string ("")
- 'cstep':       divisions of the colorbar axis ([])
- 'cmap':        colormap ("RdBu_r")
- 'proj':        projection. One of ["platecarree", "robinson", "mollweide"]. Defaults to "platecarree".
- 'cpad':        padding (shift) of the colorbar (0.08)
- 'sub':         matplotlib subplot option (e.g. "221" for the first panel of 2x2 subplots) ("111")
- 'clabel':      label of the colorbar (defaults to the units string read from the netcdf file)
- 'cdir':        direction of the colorbar. One of ["horizontal", "vertical"]. ("horizontal")
- 'cscale':      scaling of the colorbar (0.65)
- 'cfs':         colorbar ticks font size (12)
- 'lfs':         colorbar label font size (12)
- 'tfs':         title font size (14)
- 'tpad':        padding (shift) of the title string
- 'tweight':     weight of title font. One of [ 'normal" "bold" "heavy" "light" "ultrabold" "ultralight"]. ("normal")
- 'grid':        grid spacing (defaults to [60,30]). set to empty [] to remove gridlines.
- 'region::NTuple{4,Int64}':       region to plot in format (lon1, lon2, lat1, lat2). Defaults to global.
- 'style':       one of ["pcolormesh" "contourf"]. Defaults to "pcolormesh".
- 'levels':      contour plot levels. Can be an array or a number of levels (auto)
- 'extend':      plot levels outside `levels` range. One of ["neither", "both", "min", "max"]. Default: "neither".

Author: Jost von Hardenberg, 2019
"""
function plotmap(lon, lat, data; title="", cstep=[], cmap="RdBu_r", proj="",
                 cpad=0.08, tpad=24, sub=111, clabel="", cdir="horizontal",
                 cscale=0.65, tfs=14, cfs=12, lfs=12, tweight="normal",
                 grid=[60,30], region=(), style="pcolormesh", levels=0,
                 extend="neither")

dd = size(data)

if length(dd)==3 data=data[:,:,1] end
if style=="pcolormesh"
    if length(lon) in dd
        #println("pcolormesh needs cell boundaries, reconstructing lon")
        lonb=zeros(2,length(lon))
        lonb[1,2:end]=0.5*(lon[2:end]+lon[1:(end-1)])
        lonb[1,1]=lon[1]-(lon[2]-lon[1])*0.5
        lonb[2,end]=lon[end]+(lon[end]-lon[end-1])*0.5
        lon=vcat(lonb[1,:],lonb[2,end])
    end
    if length(lat) in dd
        #println("pcolormesh needs cell boundaries, reconstructing lat")
        latb=zeros(2,length(lat))
        latb[1,2:end]=0.5*(lat[2:end]+lat[1:(end-1)])
        latb[1,1]=lat[1]-(lat[2]-lat[1])*0.5
        latb[2,end]=lat[end]+(lat[end]-lat[end-1])*0.5
        if latb[1,1]>89; latb[1,1]=90 ; end
        if latb[1,1]<-89; latb[1,1]=-90 ; end
        if latb[2,end]>89; latb[2,end]=90 ; end
        if latb[2,end]<-89; latb[2,end]=-90 ; end
        lat=vcat(latb[1,:],latb[2,end])
    end
    if length(lon)==(dd[1]+1) data=data' end
else
    if length(lon)==dd[1] data=data' end
end

ccrs = pyimport("cartopy.crs")
cutil = pyimport("cartopy.util")

if proj=="robinson"
    proj=ccrs.Robinson()
    dlabels=false
elseif proj == "mollweide"
    proj=ccrs.Mollweide()
    dlabels=false
else
    proj=ccrs.PlateCarree()
    dlabels=true
end

ax = subplot(sub, projection=proj)
if length(region)>0 ax.set_extent(region, crs=ccrs.PlateCarree()) end
ax.coastlines()
xlocvec=vcat(-vcat(grid[1]:grid[1]:180)[end:-1:1], vcat(0:grid[1]:180))
ylocvec=vcat(-vcat(grid[2]:grid[2]:90)[end:-1:1], vcat(0:grid[2]:90))

if dlabels
    ax.gridlines(linewidth=1, color="gray", alpha=0.5, linestyle="--",
                 draw_labels=true, xlocs=xlocvec, ylocs=ylocvec)
else
    ax.gridlines(linewidth=1, color="gray", alpha=0.5, linestyle="--",
                 xlocs=xlocvec, ylocs=ylocvec)
end

if style=="contourf"
    data_cyc, lon_cyc = cutil.add_cyclic_point(data, coord=lon)
    if levels==0 
        contourf(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(),
                 cmap=cmap, extend=extend)
    else
        contourf(lon_cyc, lat, data_cyc, transform=ccrs.PlateCarree(),
                 cmap=cmap, levels=levels, extend=extend)
    end
else
    pcolormesh(lon, lat, data, transform=ccrs.PlateCarree(), cmap=cmap)
end

if length(cstep)>0 clim(cstep[1],cstep[end]); end
if length(title)>1 PyPlot.title(title, pad=tpad, fontsize=tfs, weight=tweight) end
cbar=colorbar(orientation=cdir, extend="both", pad=cpad, label=clabel,
              shrink=cscale)
cbar.set_label(label=clabel,size=lfs)
cbar.ax.tick_params(labelsize=cfs) 
if length(cstep)>0 cbar.set_ticks(cstep) end
tight_layout()

end
                  
