function generate_shares(δ::Matrix{Float64}, σ_α::Float64 = σ_α; p::Matrix{Float64} = p, ν_vec::Vector{Float64} = ν_vec)
    s = zeros(Float64, size(p)[1], size(p)[2])

    f_δ = zeros(Float64, size(p)[1], size(p)[2])
    for ν in ν_vec 
        for m in 1:size(p)[2] 
            pView = @view p[:,m]
            δView = @view δ[:,m]
            summation = sum(exp.(δView .- σ_α .* pView .*ν))
            f_δ[:,m] .= exp.(δView .- σ_α .* pView.*ν)./(1+ summation )
        end


        s .+= f_δ
    end
    return s/length(ν_vec)
end

function constructInstruments(X1 = X1, X2 = X2, X3 = X3, p = p, Z = Z, nJ = nJ, nM = nM)
    #Construct instruments
    #Demand instruments
    charInstrument = Array{Float64}(undef, 3, nJ, nM) 

    i = 1
    for characteristic in [X1, X2, X3]
        for m in 1:size(characteristic)[2]
            charInstrument[i,:,m] .= sum(characteristic[:,m]) .- characteristic[:,m]
        end
        i +=1
    end
    charInstrument2 = charInstrument[2,:,:]
    charInstrument3 = charInstrument[3,:,:]

    regInstruments= [vec(charInstrument2) ;;vec(charInstrument3) ]
    gmmInstruments = [vec(charInstrument2) ;;vec(charInstrument3) ;; vec(Z) ]
    xReg = [vec(X1);; vec(X2) ;; vec(X3) ;; vec(p) ]

    data = DataFrame(
        δ = vec(X1),
        X2 = vec(X2),
        X3 = vec(X3),
        p = vec(p),
        instX2 = vec(charInstrument2),
        instX3 = vec(charInstrument3)
    )
 
    return regInstruments, gmmInstruments, xReg, data

end

function innerLoop(σ_α; s_sim = s_sim, f_δ = f_δ, δ_sim = δ_sim, gmmInstruments = gmmInstruments,
        xReg = xReg, regInstruments = regInstruments, s = s, ν_vec = ν_vec, data = data)
    σ_α = σ_α[1]
    #Calculate the δ that solves the contraction mapping
    δ_sim = fixedpoint(δ -> δ .- log.(s) .+ log.(generate_shares(δ, σ_α)), zeros(Float64, nJ, nM)).zero
    data.δ = vec(δ_sim)
    ξ_sim = fit(EconometricModel,
            @formula(δ ~ X2 + X3 + (p ~ instX2 + instX3)),
            data)
    display(ξ_sim)
    ξ_sim = residuals(ξ_sim)
    #weight = inv(gmmInstruments'* ξ_sim *ξ_sim' * gmmInstruments)
    g_θ = ξ_sim' * gmmInstruments * gmmInstruments' * ξ_sim    
    return g_θ'*g_θ 



end
