__precompile__(true)

module ConvexPolicyOptimizer

using Distributed, Printf, LinearAlgebra, Statistics, Random, SparseArrays
using Parameters, ProgressMeter, Iconv, SortingAlgorithms 
using ScikitLearnBase, JuMP, OSQP, AmplNLWriter
using OSQP.OSQPMathProgBaseInterface: OSQPSolver
import ScikitLearnBase: BaseEstimator, fit!, fit_transform!, predict, transform
using MPI: Barrier, SUM, Comm_size, Comm_rank, COMM_WORLD, Initialized, Allreduce, Allgather
import MLSuiteBase: paramgrid

export ConvexOpt, CPO

const CPO = ConvexPolicyOptimizer

include("util.jl")
include("osqp.jl")
include("mpi.jl")
include("discretizer.jl")
include("admm.jl")
include("opt.jl")

function __init__()
    include(joinpath(@__DIR__, "../deps/env.jl"))
    @eval using Mosek, Gurobi
end

end