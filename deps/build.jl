using Pkg, BinDeps
using BinDeps: generate_steps, getallproviders, lower, PackageManager

!Sys.islinux() && exit()

if isnothing(Sys.which("sudo")) # in docker
    try run(`apt update`) catch end
    try run(`yum update`) catch end
end

@BinDeps.setup

wget = library_dependency("wget")
svn = library_dependency("svn")
make = library_dependency("make")
cmake = library_dependency("cmake")
gcc = library_dependency("gcc")
gfortran = library_dependency("gfortran")

common = Dict("wget" => wget, "subversion" => svn, "make" => make, "cmake" => cmake)
provides(AptGet, Dict(common..., "g++" => gcc, "gfortran" => gfortran))
provides(Yum, Dict(common..., "gcc-c++" => gcc, "gcc-gfortran" => gfortran))

for dep in bindeps_context.deps
    dp, opts = getallproviders(dep, PackageManager)[1]
    cmd = lower(generate_steps(dep, dp, opts)).steps[1]
    i = findfirst(x -> x == "install", cmd.exec)
    insert!(cmd.exec, i + 1, "-y")
    run(cmd)
end

buildsh = joinpath(@__DIR__, "build.sh")
ENV["JULIA_DEPOT_PATH"] = DEPOT_PATH[1]
run(`bash $buildsh`)

include("env.jl")
Pkg.build("Gurobi")
Pkg.build("SCIP")