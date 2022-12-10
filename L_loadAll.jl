using Parameters
using DifferentialEquations
using QuadGK
using PyPlot
using Suppressor
using JLD2
using ForwardDiff
using LinearAlgebra
using DelimitedFiles
using Distributions
using StaticArrays
using Interpolations

using JuMP
using NLopt
using Optim

# optional
using BenchmarkTools

include("L_parameters.jl")
include("L_diffODE.jl")
include("L_optimization.jl")
include("L_optimByHand.jl")
include("L_plot.jl")
