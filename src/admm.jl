using Base.Iterators: repeated

export admm_consensus

function admm_consensus(opt, dim; epochs = 100, ρ = 1.0, αr = 1.0, λ = 0.0,
                    epsabs = 1e-4, epsrel = 1e-2, μ = 10.0, cb = identity)
    # mean of x, xᵢ = z
    N, z, zpre = worldsize(), zeros(dim), zeros(dim)
    ρ = N == 1 ? 0.0 : ρ
    # Lᵨ(x, z, u) = ∑ᵢ(f(xᵢ) + ρ/2 ‖xᵢ - z + uᵢ‖₂²) + λ/2 ‖z‖₂²
    x, x̂, u = zeros(dim), zeros(dim), zeros(dim)
    for t in 1:epochs
        Initialized() && Barrier(COMM_WORLD)
        @master println("----------------------------------------")
        # primal update for xᵢ
        # xᵢ := argmin(f(xᵢ) + ρ/2 ‖xᵢ - z + uᵢ‖₂²)
        x, y = opt(z, u, ρ)
        sleep(myrank() / worldsize() / 10)
        @printf("rank: %d, obj: %.2e\n", myrank(), y)
        ȳ = allmean(y)
        @master @printf("avgobj: %.2e\n", ȳ)
        if N == 1
            cb(x)
            return x
        end
        # primal update for z with relaxation
        # z = 1/N ∑ᵢ(xᵢ + uᵢ)
        copyto!(zpre, z); fill!(z, 0)
        x̂ .= αr .* x .+ (1 - αr) .* zpre
        β = N * ρ / (λ + N * ρ)
        z .+= β .* allmean(x̂ .+ u)
        # primal residual and dual residual
        # r²: ∑ᵢ‖xᵢ - z‖₂², s²:  ρ ∑ᵢ‖z - zpre‖₂²
        r = sqrt(allsum(norm(x .- z)^2))
        s = ρ * √N * norm(z .- zpre)
        ϵr = √dim * epsabs + epsrel * max(sqrt(allsum(norm(x)^2)), √N * norm(z))
        ϵs = √dim * epsabs + epsrel * ρ * sqrt(allsum(norm(u)^2))
        @master @printf("t: %d, ρ: %.2g, r: %.2g, ϵr: %.2g, s: %.2g, ϵs: %.2g\n", t, ρ, r, ϵr, s, ϵs)
        r < ϵr && s < ϵs && return z
        # dual acent for uᵢ
        # uᵢ += xᵢ - z
        u .+= x̂ .- z
        # varying penalty parameter ρ
        τ = ifelse(r > μ * s, 2.0, ifelse(s > μ * r, 0.5, 1.0))
        ρ *= τ
        u ./= τ
        cb(z)
    end
    return z
end