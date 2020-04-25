#!/usr/bin/env julia
@info "Installing Julia dependencies"

if VERSION >= v"0.7.0-DEV.2005"
    using Pkg
end

ENV["PYTHON"] = string(ENV["CONDA_PREFIX"], "/bin/python")

Pkg.activate(@__DIR__)
Pkg.instantiate()

@info "Installed Julia packages:"
Pkg.status()
