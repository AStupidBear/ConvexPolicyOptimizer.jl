
using OSQP.OSQPMathProgBaseInterface: OSQPMathProgModel, update!, resetproblem
import OSQP.OSQPMathProgBaseInterface.MathProgBase: setquadobj!

function setquadobj!(model::OSQPMathProgModel, rowidx::Vector, colidx::Vector, quadval::Vector)
    nterms = length(quadval)
    @boundscheck length(rowidx) == nterms || error()
    @boundscheck length(colidx) == nterms || error()

    # Check if only the values have changed
    Pi, Pj, Px = findnz(model.P)
    if (rowidx == Pi) & (colidx == Pj) & !model.perform_setup
        if model.sense == :Max
            # Update matrix in MathProgBase model
            copyto!(model.P.nzval, -Px)
        else
            # Update matrix in MathProgBase model
            copyto!(model.P.nzval, Px)
        end
        # Update only nonzeros of P
        update!(model.inner, Px=model.P.nzval)

        # Reset solution status
        resetproblem(model)

        return model
    end

    # Create sparse matrix from indices
    I = vcat(Pi, rowidx, colidx)
    J = vcat(Pj, colidx, rowidx)
    V = vcat(Px, quadval, quadval)
    model.P = sparse(I, J, V, size(model.P)..., (x, y) -> y)

    # Change sign if maximizing
    if model.sense == :Max
        model.P .= .-model.P
    end

    # Need new setup when we set a new P
    model.perform_setup = true

    # Reset problem solution
    resetproblem(model)

    model
end
