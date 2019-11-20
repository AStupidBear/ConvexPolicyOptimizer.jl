using Random
using Statistics
using Convex
using SCS
using Optim
using ConvexPolicyOptimizer
using MPI

Random.seed!(1000)

!isinteractive() && MPI.Init()

n, m = 2, 200
N, M = m ÷ 2, m ÷ 2
N1, N2 = Int(0.6N), Int(0.4N)
M1, M2 = Int(0.6M), Int(0.4M)

# positive examples
Y = [1.5 .+ 0.9 .* randn(1, N1)  1.5 .+ 0.7 .* randn(1, N2);
    2 .* (randn(1, N1) .+ 1)  2 .* (randn(1, N2) .- 1)]

# negative examples
X = [-1.5 .+ 0.9 .* randn(1, M1)  -1.5 .+ 0.7 .* randn(1, M2);
    2 .* randn(1, M1) .- 2  2 .* randn(1, M2) .+ 2]

x = [X Y]
y = [ones(1, N) -ones(1, M)]
A = permutedims([y .* x; y])

p = zeros(Int, m)
p[vec(y) .== 1] = sort(rand(1:10, sum(y .== 1)))
p[vec(y) .== -1] = sort(rand(11:20, sum(y .== -1)))

# consensus
using ConvexPolicyOptimizer: myrank, worldsize, @redirect
admm_consensus(size(A, 2), λ = 1.0) do z, u, ρ
    if MPI.Initialized()
        i = myrank() + 1
        Ai = A[p .== i, :]
    else
        Ai = A
    end
    x0 = zeros(size(A, 2))
    xᵒ = if @isdefined(Convex)
        x = Convex.Variable(length(x0))
        problem = minimize(sum(max(1 - Ai * x, 0.0)) + ρ / 2 * sumsquares(x - z + u))
        @redirect devnull Convex.solve!(problem, SCSSolver())
        vec(x.value)
    else
        f = x -> sum(max.(1 .- Ai * x, 0.0)) + ρ / 2 * sum(abs2, x .- z .+ u)
        df = Optim.OnceDifferentiable(f, x0, autodiff = :forward)
        res = Optim.optimize(df, x0, BFGS())
        Optim.minimizer(res)
    end
    if MPI.Initialized()
        Ai = A[p .== (myrank() + 1), :]
        yᵒ = mean(max.(1 .- Ai * xᵒ, 0.0))
    else
        Ai = A
    end
    z = MPI.Initialized() ? z : xᵒ
    yᵒ = mean(max.(1 .- Ai * z, 0.0))
    return xᵒ, yᵒ
end

!isinteractive() && MPI.Finalize()

# julia --project -e 'using Convex, SCS, Optim, ConvexPolicyOptimizer' && \
# mpirun -np 20 julia --project test/svm.jl