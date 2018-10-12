#!/usr/bin/env julia

scriptFile=@__FILE__
scriptDir=@__DIR__
#println("file: ", scriptFile)
#println("directory: ", scriptDir)


using Pkg

println("Installing the packages.")
pkgName=in
open(scriptDir * "/julia_requirements.txt") do f
    for i in enumerate(eachline(f))

      pkgId=i[1]
      pkgName=i[2]
      println(pkgId, ": ", pkgName)
      Pkg.add(pkgName)

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
