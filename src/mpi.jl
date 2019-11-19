function allsum(x)
    if Initialized()
        x = Allreduce(x, SUM, COMM_WORLD)
    end
    return x
end

function allmean(x)
    if Initialized()
        x = Allreduce(x, SUM, COMM_WORLD)
        x = x / Comm_size(COMM_WORLD)
    end
    return x
end

function part(x::AbstractArray{T, N}, dim = -2) where {T, N}
    !Initialized() && return x
    dim = dim > 0 ? dim : N + dim + 1
    dsize = size(x, dim)
    rank = myrank()
    wsize = worldsize()
    @assert wsize <= dsize
    chunk = ceil(Int, dsize / wsize)
    is = (rank * chunk + 1):min(dsize, (rank + 1) * chunk)  
    view(x, ntuple(x -> x == dim ? is : (:), N)...)
end

part(x::Number) = x

macro master(ex)
    esc(:(myrank() == 0 && $ex))
end

myrank() = Initialized() ? Comm_rank(COMM_WORLD) : 0

worldsize() = Initialized() ? Comm_size(COMM_WORLD) : 1