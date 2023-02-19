using Random, Distributions

using BenchmarkTools, Profile
using NLsolve, Optim, LinearAlgebra
using StaticArrays

nPop = 1000
nJ = 3
nM = 100
Random.seed!(111)
include("dataGen.jl")
include("demandAnalysis.jl")


#Generate observables X, p, s, W, J
X1, X2, X3 , p , s , W , Z = generate_Data()
display(p)
display(s)

#Set initial guesses
σ_α = 1
β = @MArray [1,1,1]
α = 1
γ = @MArray [2,1,1]
ν = SArray{Tuple{nPop}}(rand(LogNormal(0,1), nPop))

#Get δ calculator up and running
δ_sim = fixedpoint(δ -> δ .- log.(s) .+ log.(generate_shares(δ, ν_vec = ν)), zeros(Float64, nJ, nM)).zero

stophere
outerRuns = 1

for outerRun in 1:outerRuns
    for innerRun in 1:innerRun
        δ 
    end
end

