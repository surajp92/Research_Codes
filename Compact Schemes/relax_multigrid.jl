include("jacobi_solver.jl")


function relax_multigrid(nx, ny, dx, dy, source, u_numerical, lambda, V)
    if flag_solver == 1
        jacobi_solver_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)
    else
        println(1)
    end
end
