function *ᶜ(A, B) # tensor contraction
    Ar = reshape(A, :, size(A, ndims(A)))
    Br = reshape(B, size(B, 1), :)
    dims = (size(A)[1:end-1]..., size(B)[2:end]...)
    C = reshape(Ar * Br, dims...)
end

function hasnan(x)
    for i in eachindex(x)
        isnan(x[i]) && return true
    end
    false
end

function part(x::AbstractArray{T, N}, dim = -2) where {T, N}
    nprocs() == 1 && return x
    dim = dim > 0 ? dim : N + dim + 1
    dsize = size(x, dim)
    csize = nworkers()
    rank = myid() - 1
    @assert csize <= dsize
    chunk = ceil(Int, dsize / csize)
    is = (rank * chunk + 1):min(dsize, (rank + 1) * chunk)  
    view(x, ntuple(x -> x == dim ? is : (:), N)...)
end

part(x::Number) = x

macro redirect(src, ex)
    src = src == :devnull ? "/dev/null" : src
    quote
        io = open($(esc(src)), "a")
        o, e = stdout, stderr
        redirect_stdout(io)
        redirect_stderr(io)
        try
            $(esc(ex)); sleep(0.01)
        finally
            flush(io); close(io)
            redirect_stdout(o)
            redirect_stderr(e)
        end
    end
end