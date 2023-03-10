




function generate_δ(p; X = X, β = β, ξ = ξ, α = α)
    return X[1,:,:] .* β[1] .+ X[2,:,:] .* β[2] .+ X[3,:,:] .* β[3] .- α .* p .+ ξ
end


function generate_f_δ(p::Matrix{Float64}; δ = δ, ν = ν, σ_α = σ_α)
    #This function is the slow down, attempt to optimize as will be useful
    f_δ = zeros(Float64, size(p)[1], size(p)[2])
    for m in 1:size(p)[2]
        f_δ[:,m] = exp.(δ[:,m] .- σ_α .* p[:,m].*ν)./(1+ sum(exp.(δ[:,m] .- σ_α .* p[:,m]*ν)) )
    end
    return f_δ
end
function generate_f_δ(p::SMatrix{nJ, nM, Float64, nJ*nM}; δ = δ, ν = ν, σ_α = σ_α)
    #This function is the slow down, attempt to optimize as will be useful
    f_δ = ones(Float64, size(p)[1], size(p)[2])
    for m in 1:size(p)[2]
        f_δ[:,m] = exp.(δ[:,m] .- σ_α .* p[:,m].*ν)./(1+ sum(exp.(δ[:,m] .- σ_α .* p[:,m]*ν)) )
    end
    return f_δ
end



function firm_behavior(p::Matrix{Float64}; ν_vec= ν_vec, MC = MC, σ_α = σ_α,β = β, X=X, α = α, ξ = ξ, return_share = false)

    partial_s_jm = zeros(Float64, size(p)[1], size(p)[2])
    s= zeros(Float64, size(p)[1], size(p)[2])

    δ = generate_δ(p, X=X, α = α, β = β, ξ = ξ)
    for i in 1:nPop
        f_δ_i = generate_f_δ(p, δ = δ, ν = ν_vec[i], σ_α = σ_α)
        partial_s_jm .+= (-α .- σ_α .* ν_vec[i]) .* (f_δ_i .* (1 .-f_δ_i)) 
        s .+= f_δ_i

    end
    
    if return_share 
       return (s/nPop, MC .- s ./ partial_s_jm )

    end
    return MC .- s ./ partial_s_jm 
end

function generate_Data()
    β = [5,1,1]
    α = 1.0
    σ_α = 1.0
    γ = [2,1,1]

    X = Array{Float64}(undef, 3, nJ, nM) 
    X[1,:,:] .= ones(size(X)[2:3])
    X[2,:,:] .= rand(Uniform(0,1), size(X)[2:3])
    X[3,:,:] .= rand(Normal(0,1), size(X)[2:3])
    ξ = rand(Normal(0,1), nJ, nM)
    ν_vec = rand(LogNormal(0,1), nPop)
    W = rand(LogNormal(0,1), nJ)
    Z = rand(LogNormal(0,1), nJ, nM)
    η = rand(LogNormal(0,1), nJ,nM)

    MC = γ[1] .+ γ[2] .* repeat(W,1,nM).+ γ[3] .* Z + η

    p = fixedpoint(p -> firm_behavior(p, X = X, ν_vec = ν_vec, MC = MC, σ_α = σ_α, α = α, β =β,
        ξ =ξ), ones(nJ, nM)).zero
    display(generate_δ(p, X = X, ξ = ξ, β = β, α = α))
    δ_cool = generate_δ(p, X = X, ξ = ξ, β = β, α = α)
    display(generate_shares(δ_cool, σ_α, p = p,
        ν_vec = ν_vec))
    s = firm_behavior(p, X = X, ν_vec = ν_vec, MC = MC, σ_α = σ_α, α = α, β =β,
        ξ =ξ, return_share = true)[1]
    X1 = X[1,:,:]
    X2 = X[2,:,:]
    X3 = X[3,:,:]
    return X1, X2, X3, p, s, W, Z, δ_cool, ν_vec 
end

