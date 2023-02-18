using Random, Distributions

using BenchmarkTools, Profile
using NLsolve, Optim, LinearAlgebra

nPop = 10000
nJ = 3
nM = 100
Random.seed!(111)
include("dataGen.jl")
include("demandAnalysis.jl")


#Generate observables X, p, s, W, J
X, p, s, W, Z = generate_Data()
display(p)
display(s)
sotpher
#
δ =  fixedpoint(δ -> δ .+ log.(s.+eps(Float64)) .- log.(generate_shares(δ).+eps(Float64)), zeros(nJ, nM))
stopher
ξ_calc = δ_calc .- X[1,:,:].*β[1]  .- X[2,:,:].*β[2] .- X[3,:,:].*β[3]

ν_vec = rand(LogitNormal(0,1), nPop)


