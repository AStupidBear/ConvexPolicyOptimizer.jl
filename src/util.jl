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

function write_feaimpt(w, c; dst = "feaimpt.csv")
    length(w) != length(c) && return
    w, c = sortall(vec(w), vec(c), by = abs, rev = true)
    open(dst, "w") do io
        for i in 1:min(100, length(c))
            write(io, togbk(c[i]), ',')
            println(io, trunc(w[i], digits = 2))
        end
    end
end

function sortall(xs::AbstractArray...; kw...)
    p = sortperm(first(xs); kw...)
    map(x -> x[p], xs)
end