include("mg2.jl")
#---------------------------- Multigrid solver---------------------------------------------------
# this function selects the multigrid solver based on number of levels to be used
# need to be implemented like a recursive function (to do)
#-------------------------------------------------------------------------------
function multigrid_solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
                initial_rms, maximum_iterations, tiny, lambda, output)
    if flag_multigrid == 2
        # @enter mg2(dx, dy, nx, ny, residual, source, u_numerical, rms,
        #     initial_rms, maximum_iterations, tiny, lambda, output)

        mg2(dx, dy, nx, ny, residual, source, u_numerical, rms,
            initial_rms, maximum_iterations, tiny, lambda, output)
    else
        println(0)
    end
end
