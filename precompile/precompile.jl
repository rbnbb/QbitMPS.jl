# Use this script to create sysimage.so for faster loading
# Execute from project root directory
# $ julia precompile/precompile.jl
# Make sure PackageCompiler is installed in the general registry
# Make sure you instantiated the root and ./test project once
# Precompilation takes about 5 minutes
using PackageCompiler

# run all tests as a precompilation script
fname = "precompilation_script.jl"
file = open(fname, "w")
write(file, "using Pkg; Pkg.test()")
close(file)

PackageCompiler.create_sysimage(
    ["ITensors", "LinearAlgebra"];
    sysimage_path = "precompile/sysimage.so",
    precompile_execution_file = fname,
    project = ".",
)

rm(fname)
