#!/usr/bin/env julia

scriptDir=@__DIR__

if VERSION >= v"0.7.0-DEV.2005"
    using Pkg
end

@info "Installing the packages from" scriptDir * "/julia_requirements.txt"
pkgName=in
open(scriptDir * "/julia_requirements.txt") do f
    for i in enumerate(eachline(f))

      pkgId=i[1]
      pkgName=i[2]
      @info "Installing" pkgName
      Pkg.add(pkgName)

      @info "Testing: ", pkgName
      # load the package this needs to be called at top-level
      Expr(:toplevel, :(module ($pkgName) end))

    end
end

# Show the package list
@info "Installed Julia packages:"
Pkg.installed()
Pkg.status()
