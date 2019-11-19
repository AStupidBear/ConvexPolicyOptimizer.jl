@with_kw mutable struct ConvexOpt <: BaseEstimator
    method::String = "LP_OSQP"
    α::Float32 = 0.95
    w::Vector{Float32} = []
end

function fit!(m::ConvexOpt, x, y; dim = 0, period = 2^32, param = nothing)
    @unpack method, α = m
    @unpack r, c⁺, c⁻, λ = y
    model = getmodel(method)
    N, T = size(r)
    F = ndims(x) == 3 ? size(x, 1) : 
             dim == 0 ? maximum(x) : dim
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
    if isnothing(param)
        ρ = 0
        @expression(model, penalty, 0)
    else
        z, u, ρ = param
        @expression(model, penalty, sum((w[f] - z[f] + u[f])^2 for f in 1:F))
    end
    @expression(model, pnl, sum(r[n, t] * s[n, t] - c⁺[n, t] * Δs⁺[n, t] - c⁻[n, t] * Δs⁻[n, t] for t in 1:T for n in 1:N))
    @objective(model, Min, (ρ / 2) * penalty - (λ / N / T) * pnl)
    log = @sprintf("JuMP-%d.log", myid())
    @redirect(log, solve(model))
    obj = getobjectivevalue(model)
    w, s = getvalue(w), getvalue(s)
    @assert length(w) == F
    @pack! m = w
    @printf("objective: %.2e\n", obj)
    return m
end

function fit_admm!(m::ConvexOpt, x, y; ka...)
    dim = ndims(x) == 2 ? maximum(x) : size(x, 1) 
    cb = function (w)
        m.w = copy(w)
        pnl = test(m, x, y)
        @printf("pnl: %.2e\n", pnl)
    end
    dst = randstring() * ".bson"
    BSON.bson(dst, x = x, y = y)
    admm_consensus(dim; cb = cb, ka...) do z, u, ρ
        data = bsload(dst, mmaparrays = true)
        x, y = part(data[:x]), map(part, data[:y])
        fit!(m, x, y, dim = dim, param = (z, u, ρ))
        hasnan(m.w) && fill!(m.w, 0)
        pnl = test(m, x, y)
        pnl < -10 && fill!(m.w, 0)
        return m.w
    end
    rm(dst, force = true)
    return m
end

function predict(m::ConvexOpt, x)
    if ndims(x) == 2
        return m.w[x]
    else
        dropdims(transpose(m.w) *ᶜ x, dims = 1)
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

function getmodel(method)
    nthreads = nprocs() > 1 ? 1 : 0
    if method == "LP_OSQP"
        solver = OSQPSolver()
    elseif method == "LP_GUROBI"
        solver = GurobiSolver(Method = 2, Crossover = 0, Threads = nthreads)
    elseif method == "LP_MOSEK"
        solver = MosekSolver(MSK_IPAR_INTPNT_BASIS = 0, MSK_IPAR_NUM_THREADS = nthreads)
    elseif method == "NLP"
        ka = Dict(:linear_solver => "ma86", :max_cpu_time => 3e3)
        @static Sys.iswindows() && delete!(options, :linear_solver)
        options = [string(k) * "=" * string(v) for (k, v) in ka]
        solver = AmplNLSolver(IPOPT_AMPLEXE, options)
    else
        error("no solver loaded")
    end
    m = Model(solver = solver)
end