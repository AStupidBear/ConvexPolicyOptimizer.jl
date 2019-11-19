using ConvexPolicyOptimizer
using Test
const CPO = ConvexPolicyOptimizer

F, N, T = 10, 10, 100
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
    CPO.fit_admm!(model, x, y)
    CPO.fit_admm!(model, x_int, y)
end