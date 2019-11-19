using Random
using ConvexPolicyOptimizer
using MPI
using Test

!isinteractive() && MPI.Init()
cd(tempdir())

F, N, T = 10, 50, 100
Random.seed!(1234)
x = randn(Float32, F, N, T)
x_int = rand(0x01:0x10, N, T)
y = @. (x[1, :, :] > 0) & (x[2, :, :] > 0)
r = @. ifelse(y > 0, 1f0, -0.2f0)
c⁺ = c⁻ = fill(1f-4, N, T)
y = (r = r, c⁺ = c⁺, c⁻ = c⁻, λ = 1f0)

for method in ["LP_OSQP", "LP_MOSEK", "NLP"]
    model = ConvexOpt(method = method)
    CPO.fit!(model, x, y)
    CPO.fit!(model, x_int, y)
end

!isinteractive() && MPI.Finalize()

# julia --project -e 'using ConvexPolicyOptimizer' && \
# mpirun -np 10 julia --project test/runtests.jl