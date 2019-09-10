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
