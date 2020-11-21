module ModifiedLatinHypercubeSampling
# Hess, Stephane, Kenneth E. Train and John W. Polak, (2006), "On the Use of a Modified Latin Hypercube
# Sampling (MLHS) Method in the Estimation of a Mixed Logit Model for Vehicle Choice,"
# Transportation Research Part B, Vol. 40, pp. 147 - 163.

using Random
export mlhs

# R draws for D dimensions
# R is the number of draws
# D is the number of dimensions

function mlhs(R, D)

    # ϕ from Equation (14)
    ϕ = (0:(R-1)) ./ R
    ϕ_grid = repeat(ϕ, 1, D)

    # Obtain random value
    Random.seed!(19900729)
    randval = rand(D)
    while 0.0 in randval || 1.0 in randval
        randval = rand(D)
    end
    x = transpose(randval / R)

    # Ψ from Equation (15)
    Ψ = ϕ .+ x

    # Shuffle sequences in each dimension
    Ψ_shuffled = zeros(R, D)

    for i in 1:D
        Ψ_shuffled[:,i] = Ψ[shuffle!(collect(1:R)),i]
    end
    Ψ_shuffled
end

end # Module

