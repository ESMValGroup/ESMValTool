#!/usr/bin/env julia

scriptFile=@__FILE__
scriptDir=@__DIR__
#println("file: ", scriptFile)
#println("directory: ", scriptDir)

if VERSION >= v"0.7.0-DEV.2005"
    using Pkg
end

println("Installing the packages.")
pkgName=in
open(scriptDir * "/julia_requirements.txt") do f
    for i in enumerate(eachline(f))

      pkgId=i[1]
      pkgName=i[2]

      println(pkgId, ": ", pkgName)
if VERSION >= v"0.7.0-DEV.2005"
      if occursin("https://", pkgName)
          Pkg.clone(pkgName)
      else
          Pkg.add(pkgName)
      end
else
      if contains(pkgName, "https://")
          Pkg.clone(pkgName)
      else
          Pkg.add(pkgName)
      end
end
      println("Testing: ", pkgName)
      # load the package this needs to be called at top-level
      Expr(:toplevel, :(module ($pkgName) end))

    end
end

# Update the packages
#Pkg.update()

# Show the package list
println("Installed Julia packages:")
Pkg.installed()
Pkg.status()
