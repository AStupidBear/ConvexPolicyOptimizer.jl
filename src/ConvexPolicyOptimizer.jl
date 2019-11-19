__precompile__(true)

module ConvexPolicyOptimizer

using Distributed, Printf, LinearAlgebra, Statistics, Random 
using Parameters, BSON, BSONMmap, ScikitLearnBase 
using JuMP, OSQP, AmplNLWriter
using OSQP.OSQPMathProgBaseInterface: OSQPSolver
import ScikitLearnBase: BaseEstimator, fit!, predict

export ConvexOpt, fit_admm!

include("util.jl")
include("osqp.jl")
include("admm.jl")
include("opt.jl")

function __init__()
    include(joinpath(@__DIR__, "../deps/env.jl"))
    @eval using Mosek, Gurobi
end

end