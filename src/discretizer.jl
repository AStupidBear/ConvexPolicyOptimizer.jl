@with_kw mutable struct KBinsDiscretizer
    n_bins::Int = 8
    bin_edges::Vector{Vector{Float32}} = []
end

function fit!(disc::KBinsDiscretizer, x)
    resize!(disc.bin_edges, size(x, 1))
    p = range(0, 1, length = disc.n_bins + 1)[2:(end - 1)]
    disc.bin_edges = @showprogress pmap(eachrow(x)) do xf
        sort!(xf, alg = RadixSort)
        quantile(xf, p, sorted = true)
    end
end

function transform(disc::KBinsDiscretizer, x)
    xd = BitArray(undef, disc.n_bins, size(x))
    @inbounds for n in 1:size(x, 2), f in 1:size(x, 1)
        edge = disc.bin_edges[f]
        i = searchsortedfirst(edge, x[f, n])
        xd[i, f, n] = 1
    end
    reshape(xd, :, size(xd)[end])
end