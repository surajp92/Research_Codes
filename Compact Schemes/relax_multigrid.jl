include("jacobi_solver.jl")
include("gauss_seidel.jl")
include("steepest_descent.jl")
include("conjugate_gradient.jl")
include("biconjugate_gradient_stab.jl")
include("jacobi_solver_compact.jl")
include("gauss_seidel_compact.jl")
include("steepest_descent_compact.jl")
include("conjugate_gradient_compact.jl")
include("biconjugate_gradient_stab_compact.jl")


#------------------------ Relaxation multigrid----------------------------------
# this function select ths relaxation step solver based on flag_solver
#
#-------------------------------------------------------------------------------
function relax_multigrid(nx, ny, dx, dy, source, u_numerical, lambda, tiny,
                         V, flag_solver, flag_order)
    if flag_order == 1
        if flag_solver == 1
    # uncomment below line if debugging using ASTInterpreter2
    #    @ enter jacobi_solver_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)
            jacobi_solver_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)
        elseif flag_solver == 2
            gauss_seidel_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)
        elseif flag_solver == 3
            steepest_descent_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)
        elseif flag_solver == 4
            conjugate_gradient_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)
        elseif flag_solver == 5
            biconjugate_gradient_stab_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)
        end

    elseif flag_order == 2
        # println("Compact scheme in progress")
        if flag_solver == 1
            jacobi_solver_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)
        # elseif flag_solver == 2
        #     gauss_seidel_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)
        # elseif flag_solver == 3
        #     steepest_descent_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)
        # elseif flag_solver == 4
        #     conjugate_gradient_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)
        # elseif flag_solver == 5
        #     biconjugate_gradient_stab_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)
        end
    end
end
