include("gauss_seidel.jl")
include("jacobi_solver.jl")
include("steepest_descent.jl")
include("conjugate_gradient.jl")
include("biconjugate_gradient_stab.jl")




function solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
                initial_rms, maximum_iterations, tiny, lambda)
#-------------------------------------------------------------------------------
# This function selects the iterative solver based on the flag_solver
# Arguments:
#   nx, ny = number of grid in x and y direction
#   dx, dy = grid spacing in x and y direction
#   source = source term on rhs of poisson equation ∇2ϕ = S
#   u_numerical = numerical solution for different problems
#   residual = residual at each grid location
#   rms, initial_rms = root mean squire of residual
#   maximum_iterations = maximum number of iterations
#   tiny = very small number to avoid division by zero
#-------------------------------------------------------------------------------
    if flag_multigrid == 1
        if flag_solver == 1
        # call jacobi solver
            jacobi_solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
                          initial_rms, maximum_iterations, lambda)
        elseif flag_solver == 2
        # call gauss seidel solver
            gauss_seidel(dx, dy, nx, ny, residual, source, u_numerical, rms,
                         initial_rms, maximum_iterations, lambda)
        elseif flag_solver == 3
        # call steepest descent solver
            steepest_descent(dx, dy, nx, ny, residual, source, u_numerical, rms,
                             initial_rms, maximum_iterations, tiny, lambda)
        elseif flag_solver == 4
        # call conjugate gradient solver
            conjugate_gradient(dx, dy, nx, ny, residual, source, u_numerical, rms,
                               initial_rms, maximum_iterations, tiny, lambda)
        else
            biconjugate_gradient_stab(dx, dy, nx, ny, residual, source, u_numerical,
                                      rms, initial_rms, maximum_iterations, tiny, lambda)
        end
    else
        multigrid_solver()
    end
end
