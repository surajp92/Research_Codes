include("mg2.jl")


function multigrid_solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
                initial_rms, maximum_iterations, tiny, lambda, output)
    if flag_multigrid == 2
        mg2(dx, dy, nx, ny, residual, source, u_numerical, rms,
            initial_rms, maximum_iterations, tiny, lambda, output)
    else
        println(0)
    end
end
