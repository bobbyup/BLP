function generate_shares(δ::Matrix{Float64}; p = p, X = X, ν_vec = ν_vec)
    s = zeros(Float64, size(p)[1], size(p)[2])
    for i in 1:length(ν_vec)
        s .+=generate_f_δ(p, δ = δ, ν = ν_vec[i])
    end
    return s/length(ν_vec)
end

