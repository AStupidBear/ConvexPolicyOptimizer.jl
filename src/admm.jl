using Base.Iterators: repeated

export admm_consensus

function admm_consensus(opt, dim; epochs = 100, njobs = nworkers(), ρ = 1.0, αr = 1.0,
            λ = 0.0, epsabs = 1e-4, epsrel = 1e-2, μ = 10.0, cb = identity)
    # mean of x, xᵢ = z
    N, z, zpre = njobs, zeros(dim), zeros(dim)
    ρ = N == 1 ? 0.0 : ρ
    # Lᵨ(x, z, u) = ∑ᵢ(f(xᵢ) + ρ/2 ‖xᵢ - z + uᵢ‖₂²) + λ/2 ‖z‖₂²
    xs = [zeros(dim) for n in 1:N]
    x̂s = [zeros(dim) for n in 1:N]
    us = [zeros(dim) for n in 1:N]
    for t in 1:epochs
        # primal update for xᵢ
        # xᵢ := argmin(f(xᵢ) + ρ/2 ‖xᵢ - z + uᵢ‖₂²)
        zs, ρs = repeated(z, N), repeated(ρ, N)
        xs .= map(opt, zs, us, ρs)
        N == 1 && return first(xs)
        # primal update for z with relaxation
        # z = 1/N ∑ᵢ(xᵢ + uᵢ)
        copyto!(zpre, z); fill!(z, 0)
        for (x, x̂) in zip(xs, x̂s)
            x̂ .= αr .* x .+ (1 - αr) .* zpre
        end
        β = N * ρ / (λ + N * ρ) / N
        for (x̂, u) in zip(x̂s, us)
            z .+= β .* (x̂ .+ u)
        end
        # primal residual and dual residual
        # r²: ∑ᵢ‖xᵢ - z‖₂², s²:  ρ ∑ᵢ‖z - zpre‖₂²
        r = sqrt(sum(norm(x .- z)^2 for x in xs))
        s = ρ * √N * norm(z .- zpre)
        ϵr = √dim * epsabs + epsrel * max(sqrt(sum(norm(x)^2 for x in xs)), √N * norm(z))
        ϵs = √dim * epsabs + epsrel * ρ * sqrt(sum(norm(u)^2 for u in us))
        @printf("t: %d, ρ: %.2g, r: %.2g, ϵr: %.2g, s: %.2g, ϵs: %.2g\n", t, ρ, r, ϵr, s, ϵs)
        r < ϵr && s < ϵs && return z
        # dual acent for uᵢ
        # uᵢ += xᵢ - z
        for (u, x̂) in zip(us, x̂s)
            u .+= x̂ .- z
        end
        # varying penalty parameter ρ
        τ = ifelse(r > μ * s, 2.0, ifelse(s > μ * r, 0.5, 1.0))
        ρ *= τ
        for u in us
            u ./= τ
        end
        cb(z)
    end
    return z
end