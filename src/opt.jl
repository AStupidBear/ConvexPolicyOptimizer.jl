@with_kw mutable struct ConvexOpt <: BaseEstimator
    w::Vector{Float32} = []
    disc::KBinsDiscretizer = KBinsDiscretizer()
    method::String = "LP_OSQP"
    α::Float32 = 0.95
    n_bins::Int = 8
end

paramgrid(::ConvexOpt) = [Dict("α" => a) for a in [0.5:0.1:0.8..., 0.9:0.01:0.99..., 2, 3]]

function _fit!(m::ConvexOpt, x, y, z, u, ρ; period = 2^32)
    @unpack method, α = m
    @unpack r, c⁺, c⁻, λ = y
    model = getmodel(method)
    F, N, T = length(z), size(r, 1), size(r, 2)
    # action/signal
    if ndims(x) == 2
        @variable(model, -1 <= w[f = 1:F] <= 1)
        if startswith(method, "NLP")
            @NLexpression(model, a[n = 1:N, t = 1:T], w[x[n, t]])
        else
            @expression(model, a[n = 1:N, t = 1:T], w[x[n, t]])
        end
    else
        @variable(model, a[n = 1:N, t = 1:T])
        @variable(model, -5 <= w[f = 1:F] <= 5)
        for t in 1:T, n in 1:N
            @constraint(model, a[n, t] == sum(w[f] * x[f, n, t] for f in 1:F))
        end
    end
    # state/postition
    @variable(model, -1 <= s[n = 1:N, t = 1:T] <= 1)
    if startswith(method, "LP")
        for t in 1:T, n in 1:N
            if t % period == 0 || t == 1
                @constraint(model, s[n, t] == a[n, t])
            else
                @constraint(model, s[n, t] == α * s[n, t - 1] + (1 - α) * a[n, t])
            end
        end
    elseif startswith(method, "NLP")
        for t in 1:T, n in 1:N
            if t % period == 0 || t == 1
                @NLconstraint(model, s[n, t] == tanh(a[n, t]))
            else
                @NLconstraint(model, s[n, t] == tanh(α * s[n, t - 1] + a[n, t]))
            end
        end
    end
    # slack variable
    @variable(model, Δs⁺[n = 1:N, t = 1:T] >= 0)     # Δs⁺[t] = max(s[t] - s[t - 1], 0)
    @variable(model, Δs⁻[n = 1:N, t = 1:T] >= 0)     # Δs⁻[t] = max(s[t - 1] - s[t], 0)
    for t in 1:T, n in 1:N
        if t == 1
            @constraint(model, s[n, t] <= Δs⁺[n, t])
            @constraint(model, s[n, t] >= -Δs⁻[n, t])
        else
            @constraint(model, s[n, t] - s[n, t - 1] <= Δs⁺[n, t])
            @constraint(model, s[n, t] - s[n, t - 1] >= -Δs⁻[n, t])
        end
    end
    # objective
    @expression(model, penalty, sum((w[f] - z[f] + u[f])^2 for f in 1:F))
    @expression(model, pnl, sum(r[n, t] * s[n, t] - c⁺[n, t] * Δs⁺[n, t] - c⁻[n, t] * Δs⁻[n, t] for t in 1:T for n in 1:N))
    @objective(model, Min, (ρ / 2) * penalty - (λ / N / T) * pnl)
    log = @sprintf("JuMP-%d.log", myrank())
    @redirect(log, solve(model))
    m.w = getvalue(w)
    pnl = getvalue(pnl)
    return (λ / N / T) * pnl
end

function fit!(m::ConvexOpt, x, y; columns = nothing, ka...)
    @unpack n_bins = m
    m.disc = KBinsDiscretizer(n_bins = n_bins)
    if ndims(x) == 2
        dim = maximum(x)
        columns = string.("leaf:", 1:dim)
    else
        columns = something(columns, string.(1:size(x, 1)))
        x = fit_transform!(m.disc, x)
        dim = size(x, 1)
    end
    cb = function (w)
        @pack! m = w
        x′, y′ = part(x), map(part, y)
        pnl = allmean(test(m, x′, y′))
        @master @printf("avgpnl: %.2e\n", pnl)
        @master visualize(m, columns)
    end
    admm_consensus(dim; cb = cb, ka...) do z, u, ρ
        x′, y′ = part(x), map(part, y)
        pnl = _fit!(m, x′, y′, z, u, ρ)
        @printf("rank: %d, pnl: %.2e\n", myrank(), pnl)
        hasnan(m.w) && fill!(m.w, 0)
        pnl = test(m, x′, y′)
        pnl < -10 && fill!(m.w, 0)
        return m.w
    end
    return m
end

function predict(m::ConvexOpt, x)
    if ndims(x) == 2
        return m.w[x]
    else
        xd = transform(m.disc, x)
        dropdims(transpose(m.w) *ᶜ xd, dims = 1)
    end
end

function test(m::ConvexOpt, x, y)
    @unpack α, method, w = m
    @unpack r, c⁺, c⁻, λ = y
    ŷ = predict(m, x)
    trans = if startswith(method, "LP")
        (s, a) -> α * s + (1 - α) * a
    elseif startswith(method, "NLP")
        (s, a) -> tanh(α * s + a)
    end
    N, T = size(r, 1), size(r, 2)
    pnl, p = 0.0, zeros(N)
    for t in 1:T, n in 1:N
        s, a = p[n], ŷ[n, t]
        s′ = clamp(trans(s, a), -1, 1)
        Δs = s′ - s
        Δs⁺ = max(Δs, 0)
        Δs⁻ = max(-Δs, 0)
        pnl += r[n, t] * s′ - c⁺[n, t] * Δs⁺ + c⁻[n, t] * Δs⁻
        p[n] = s′
    end
    return λ / N / T * pnl 
end

function visualize(m::ConvexOpt, columns)
    @unpack w, disc = m
    if length(disc.bin_edges) > 0
        columns = ["$c:$n" for c in columns for n in 1:disc.n_bins]
    end
    write_feaimpt(w, columns)
end

function getmodel(method)
    nthreads = worldsize() > 1 ? 1 : 0
    if method == "LP_OSQP"
        solver = OSQPSolver()
    elseif method == "LP_GUROBI"
        solver = GurobiSolver(Method = 2, Crossover = 0, Threads = nthreads)
    elseif method == "LP_MOSEK"
        solver = MosekSolver(MSK_IPAR_INTPNT_BASIS = 0, MSK_IPAR_NUM_THREADS = nthreads)
    elseif method == "NLP"
        worldsize() > 1 && (ENV["OMP_NUM_THREADS"] = 1)
        ka = Dict(:linear_solver => "ma86", :max_cpu_time => 3e3)
        @static Sys.iswindows() && delete!(options, :linear_solver)
        options = [string(k) * "=" * string(v) for (k, v) in ka]
        solver = AmplNLSolver(IPOPT_AMPLEXE, options)
    else
        error("no solver loaded")
    end
    m = Model(solver = solver)
end