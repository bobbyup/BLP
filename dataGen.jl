




function generate_δ(p::Matrix{Float64}; X = X, β = β, ξ = ξ, α = α)
    return X[1,:,:] .* β[1] .+ X[2,:,:] .* β[2] .+ X[3,:,:] .* β[3] .- α .* p .+ ξ
end

function generate_f_δ(p::Matrix{Float64}; δ = δ, ν = ν, σ_α = σ_α)
    #This function is the slow down, attempt to optimize as will be useful
    f_δ = ones(Float64, size(p)[1], size(p)[2])
    for m in 1:size(p)[2]
        f_δ[:,m] = exp.(δ[:,m] .- σ_α .* p[:,m].*ν)./(1+ sum(exp.(δ[:,m] .- σ_α .* p[:,m]*ν)) +eps(Float64))
    end
    return f_δ
end


function firm_behavior(p::Matrix{Float64}; ν_vec= ν_vec, MC = MC, σ_α = σ_α, X=X, return_share = false)

    partial_s_jm = zeros(Float64, size(p)[1], size(p)[2])
    s= zeros(Float64, size(p)[1], size(p)[2])

    δ = generate_δ(p, X=X)
    for i in 1:nPop
        f_δ_i = generate_f_δ(p, δ = δ, ν = ν_vec[i])
        partial_s_jm .+= (-α .- σ_α .* ν_vec[i]) .* (f_δ_i .* (1 .-f_δ_i)) 
        s .+= f_δ_i

    end
    
    partial_s_jm .+= partial_s_jm .+ eps(Float32)
    if return_share 
       return (s/nPop, MC .+ s ./ partial_s_jm )

    end
    return MC .- s ./ partial_s_jm 
end

function generate_Data()
    β = [5,1,1]
    α = 1
    σ_α = 1
    γ = [2,1,1]

    X = Array{Float64}(undef, 3, nJ, nM) 
    #Make the rows j and columns m
    X[1,:,:] .= ones(size(X)[2:3])
    X[2,:,:] .= rand(Uniform(0,1), size(X)[2:3])
    X[3,:,:] .= rand(Normal(0,1), size(X)[2:3])



    ξ = rand(Normal(0,1), nJ, nM)
    ν_vec = rand(LogitNormal(0,1), nPop)
    W = rand(LogNormal(0,1), nJ)
    Z = rand(LogNormal(0,1), nJ, nM)
    η = rand(LogitNormal(0,1), nJ,nM)

    
    MC .= γ[1] .+ γ[2] .* W .+ γ[3] .* Z + η

    p = fixedpoint(p -> firm_behavior(p, X = X), ones(nJ, nM)).zero

    s = firm_behavior(p, return_share = true)[1]
    return X, p, s, W, Z
end

