using Random, Distributions
using Econometrics
using BenchmarkTools, Profile
using NLsolve, Optim, LinearAlgebra
using StaticArrays
using DataFrames
using StatsBase
using Optim
using Tullio
nPop = 100
nJ = 3
nM = 100
Random.seed!(111)
include("dataGen.jl")
include("demandAnalysis.jl")


#Generate observables X, p, s, W, J
X1, X2, X3 , p , s , W , Z , δ_cool, ν_vec = generate_Data()

#@btime generate_shares(δ_cool, 1.0)
display(s)
#Set initial guesses
σ_α = 1.0
α = 1.0

ν_vec = rand(LogNormal(0,1), nPop)


#Get δ calculator up and running
s_sim = zeros(Float64, size(p)[1], size(p)[2])
f_δ = zeros(Float64, size(p)[1], size(p)[2])
δ_sim =zeros(Float64, size(p)[1], size(p)[2])

 
#Construct Instruments and regression data set
regInstruments , gmmInstruments, xReg, data = constructInstruments()
δ_sim = innerLoop([1.0])
