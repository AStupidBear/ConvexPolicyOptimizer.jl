@with_kw mutable struct KBinsDiscretizer <: BaseEstimator
    n_bins::Int = 8
    bin_edges::Vector{Vector{Float32}} = []
end

function fit!(disc::KBinsDiscretizer, x, y = nothing)
    eltype(x) <: Integer && return disc
    x = reshape(x, size(x, 1), :)
    resize!(disc.bin_edges, size(x, 1))
    p = range(0, 1, length = disc.n_bins + 1)[2:(end - 1)]
    pp = range(0, 1, length = 10 * disc.n_bins + 1)[2:(end - 1)]
    prog = Progress(size(x, 1), desc = "fitting discretizer...")
    for f in 1:size(x, 1)
        xf = x[f, :]
        sort!(xf, alg = RadixSort)
        if Initialized()
            v = quantile(xf, pp, sorted = true)
            v = Allgather(v, COMM_WORLD)
            disc.bin_edges[f] = quantile(v, p)
        else
            disc.bin_edges[f] = quantile(xf, p, sorted = true)
        end
        @master next!(prog)
    end
    return disc
end

function transform(disc::KBinsDiscretizer, x)
    eltype(x) <: Integer && return x
    sizes = size(x)
    x = reshape(x, size(x, 1), :)
    xd = BitArray(undef, disc.n_bins, size(x)...)
    @inbounds for n in 1:size(x, 2), f in 1:size(x, 1)
        edge = disc.bin_edges[f]
        i = searchsortedfirst(edge, x[f, n])
        xd[i, f, n] = true
    end
    reshape(xd, size(xd, 1) * sizes[1], sizes[2:end]...)
end