using Base.Iterators: repeated
using JuMP, OSQP, Ipopt, AmplNLWriter
using OSQP.OSQPMathProgBaseInterface: OSQPSolver

export fitadmm

function getmodel(method)
    nthreads = nprocs() > 1 ? 1 : 0
    if isdefined(Main, :Gurobi) && (method == :LP || method == :MIP)
        solver = Main.GurobiSolver(Method = 2, Crossover = 0, Threads = nthreads)
    elseif isdefined(Main, :Mosek) && method == :LP
        solver = Main.MosekSolver(MSK_IPAR_INTPNT_BASIS = 0, MSK_IPAR_NUM_THREADS = nthreads)
    elseif isdefined(Main, :OSQP) && method == :LP
        solver = Main.OSQPMathProgBaseInterface.OSQPSolver()
    elseif isdefined(Main, :AmplNLWriter) && method == :NLP
        kws = Dict(:linear_solver => "ma86", :max_cpu_time => 3e3)
        @static Sys.iswindows() && delete!(options, :linear_solver)
        options = [string(k) * "=" * string(v) for (k, v) in kws]
        solver = Main.AmplNLSolver(Ipopt.amplexe, options)
    elseif isdefined(Main, :SCIP) && method == :MIP
        solver = Main.SCIPSolver("lp/initalgorithm", 'b')
    else
        error("no solver loaded")
    end
    m = Model(solver = solver)
end

function solvemodel(m)
    src = @sprintf("JuMP-%d.log", myid())
    nprocs() > 1 ? @redirect(src, solve(m)) : solve(m)
    @printf("objective: %.2e\n", getobjectivevalue(m))
end

function fitlp(x, r, c⁺, c⁻, z, u, ρ; trunc = typemax(Int), α = 0.95, method = :LP)
    m, N, T, F = getmodel(method), size(r, 1), size(r, 2), length(z)
    # score
    if size(x, 1) == 1
        @variable(m, -1 <= w[f = 1:F] <= 1)
        @expression(m, penalty, ρ / 2 * sum((w[f] - z[f] + u[f])^2 for f in 1:F))
        if method == :NLP
            @NLexpression(m, s[n = 1:N, t = 1:T], w[Int(x[1, n, t])])
        else
            @expression(m, s[n = 1:N, t = 1:T], w[Int(x[1, n, t])])
        end
    else
        @variable(m, s[n = 1:N, t = 1:T])
        @variable(m, -5 <= w[f = 1:F] <= 5)
        @expression(m, penalty, ρ / 2 * sum((w[f] - z[f] + u[f])^2 for f in 1:F))
        for t in 1:T, n in 1:N
            @constraint(m, s[n, t] == sum(w[f] * x[f, n, t] for f in 1:F))
        end
    end
    # position
    @variable(m, p[n = 1:N, t = 1:T])
    if method == :LP
        # position
        for t in 1:T, n in 1:N
            if t % trunc == 0 || t == 1
                @constraint(m, p[n, t] == s[n, t])
            else
                @constraint(m, p[n, t] == (1 - α) * s[n, t] + α * p[n, t - 1])
            end
        end
    elseif method == :NLP
        for t in 1:T, n in 1:N
            if t % trunc == 0 || t == 1
                @NLconstraint(m, p[n, t] == tanh(s[n, t]))
            else
                @NLconstraint(m, p[n, t] == tanh(s[n, t] + α * p[n, t - 1]))
            end
        end
    elseif method == :MIP && T <= 3000 && N >= 500
        @expression(m, penalty, 0)
        @variable(m, z[ni = 1:N, nj = 1:N, t = 1:T], Bin)
        @variable(m, p̃[n = 1:N, t = 1:T], Bin)
        for t in 1:T, nj in 1:N, ni in 1:N
            @constraint(m, z[ni, nj, t] <= 1 + s[ni, t] - s[nj, t])
        end
        for t in 1:T, n in 1:N
            if t % trunc == 0 || t == 1
                @constraint(m, p[n, t] == p̃[n, t])
            else
                @constraint(m, p[n, t] == (1 - α) * p̃[n, t]  + α * p[n, t - 1])
            end
            @constraint(m, p̃[n, t] <= sum(z[n, nj, t] for nj in 1:N) / 0.8N)
        end
    elseif method == :MIP && size(x, 1) == 1
        @variable(m, z⁰[1:F], Bin)
        @variable(m, z⁺[1:F], Bin)
        for t in 1:T, n in 1:N
            f = Int(x[1, n, t])
            if t % trunc == 0 || t == 1
                @constraint(m, p[n, t] == z⁰[f])
            else
                @constraint(m, p[n, t] == p[t - 1] * z⁺[f] + (1 - p[t - 1]) * z⁰[f])
            end
        end
    elseif method == :MIP
        @variable(m, z⁻[n = 1:N, t = 1:T], Bin)
        @variable(m, z⁰[n = 1:N, t = 1:T], Bin)
        @variable(m, z⁺[n = 1:N, t = 1:T], Bin)
        @variable(m, p⁻[n = 1:N, t = 1:T])
        @variable(m, p⁰[n = 1:N, t = 1:T])
        @variable(m, p⁺[n = 1:N, t = 1:T])
        for t in 1:T, n in 1:N
            @constraint(m, z⁻[n, t] + z⁰[n, t] + z⁺[n, t] == 1)
            @constraint(m, p⁻[n, t] <= -z⁻[n, t])
            @constraint(m, p⁻[n, t] >= -5z⁻[n, t])
            @constraint(m, p⁰[n, t] <= z⁰[n, t])
            @constraint(m, p⁰[n, t] >= -z⁰[n, t])
            @constraint(m, p⁺[n, t] <= 5z⁺[n, t])
            @constraint(m, p⁺[n, t] >= z⁺[n, t])
            if t % trunc == 0 || t == 1
                @constraint(m, p⁻[n, t] + p⁰[n, t] + p⁺[n, t] == s[n, t])
            else
                @constraint(m, p⁻[n, t] + p⁰[n, t] + p⁺[n, t] == s[n, t] + α * p[n, t - 1])
            end
            @constraint(m, p[n, t] == -z⁻[n, t] + z⁺[n, t])
        end
    end
    # slack position
    @variable(m, Δp⁺[n = 1:N, t = 1:T] >= 0)     # Δp⁺[t] = max(p[t] - p[t - 1], 0)
    @variable(m, Δp⁻[n = 1:N, t = 1:T] >= 0)     # Δp⁻[t] = max(p[t] - p[t - 1], 0)
    for t in 1:T, n in 1:N
        if t == 1
            @constraint(m, p[n, t] <= Δp⁺[n, t])
            @constraint(m, p[n, t] >= -Δp⁻[n, t])
        else
            @constraint(m, p[n, t] - p[n, t - 1] <= Δp⁺[n, t])
            @constraint(m, p[n, t] - p[n, t - 1] >= -Δp⁻[n, t])
        end
    end
    #  loss
    if method == :LP && eltype(x) != Int
        @variable(m, l[n = 1:N, t = 1:T])            # loss l[t] = max(-r[t] * p[t], -abs(r[t]))
        for t in 1:T, n in 1:N
            @constraint(m, l[n, t] >= -1)
            @constraint(m, l[n, t] >= -sign(r[n, t]) * p[n, t])
        end
        @objective(m, Min, penalty + sum(abs(r[n, t]) * l[n, t] + c⁺[n, t] * Δp⁺[n, t] + c⁻[n, t] * Δp⁻[n, t] for t in 1:T for n in 1:N))
    else
        @objective(m, Min, penalty + sum(-r[n, t] * p[n, t] + c⁺[n, t] * Δp⁺[n, t] + c⁻[n, t] * Δp⁻[n, t] for t in 1:T for n in 1:N))
    end
    solvemodel(m)
    wv, pv = getvalue(w), getvalue(p)
end

function fitadmm(dim; kws...)
    cb = w -> @printf("avgloss: %.2e\n", testlp(w; kws...))
    admm_consensus(partfitadmm, dim; cb = cb, kws...)
end

function partfitadmm(z, u, ρ, i; kws...)
    x, r, c⁺, c⁻ = getlpdata("data_train.h5")
    w, p = fitlp(x, r, c⁺, c⁻, z, u, ρ; kws...)
    if hasnan(w)
        @warn("w has nans, filling w with 0")
        fill!(w, 0)
    end
    l = testlp(w; kws...)
    if l > 10
        @warn("loss > 10, filling w with 0")
        fill!(w, 0)
    end
    return w
end

function getlpdata(src)
    data = loaddata(src)
    β = 224f0 / ndays(data) / nstocks(data)
    x = part(data.特征)
    β *= 100 * nticks(data) / size(x, 3)
    r = β * part(data.涨幅)
    c⁺ = β * part(data.买手续费率)
    c⁻ = β * part(data.卖手续费率)
    return x, r, c⁺, c⁻
end

function testlp(w; kws...)
    x, r, c⁺, c⁻ = getlpdata("data_train.h5")
    testlp(w, x, r, c⁺, c⁻; kws...)
end

function testlp(w, x, r, c⁺, c⁻; α = 0.95, method = :LP, kws...)
    N, T = size(r, 1), size(r, 2)
    l, pold = 0.0, zeros(N)
    for t in 1:T, n in 1:N
        s = size(x, 1) == 1 ? w[Int(x[1, n, t])] : dot(w, x[:, n, t])
        if method == :LP
            p = (1f0 - α) * s + α * pold[n]
        elseif method == :MIP
            p = clamp(s + α * pold[n], -1, 1)
        elseif method == :NLP
            p = tanh(s + α * pold[n])
        end
        pc = clamp(p, -1, 1)
        Δp = pc - pold[n]
        Δp⁺ = max(Δp, 0)
        Δp⁻ = max(-Δp, 0)
        l += -r[n, t] * pc + c⁺[n, t] * Δp⁺ + c⁻[n, t] * Δp⁻
        pold[n] = p
    end
    return l
end
