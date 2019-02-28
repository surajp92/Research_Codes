include("residualcalculation.jl")
include("relax_multigrid.jl")
include("mg_operation.jl")

#-------------------------Multigrid Solver with 2 levels------------------------
# this function executes multigrid framework with two levels
#
#-------------------------------------------------------------------------------
function mg2(dx, dy, nx, ny, residual, source, u_numerical, rms,
            initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0
    # compute initial residual
    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)
    # compute initial L-2 norm
    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    println("0", " ", rms, " ", rms/initial_rms)
<<<<<<< HEAD:Compact Schemes/extra files/mg2.jl

=======
>>>>>>> 57e40ec981309fb68a0525bf39a226975e7b9c25:Compact Schemes/mg2.jl

    #allocate memory for grid size at different levels
    level_nx = zeros(Int64, 2)
    level_ny = zeros(Int64, 2)
    level_dx = zeros(Float64, 2)
    level_dy = zeros(Float64, 2)
    level_nx[1] = nx
    level_ny[1] = ny
    level_nx[2] = Int64(level_nx[1]/2)
    level_ny[2] = Int64(level_ny[1]/2)
    level_dx[1] = dx
    level_dy[1] = dy
    level_dx[2] = level_dx[1]*2
    level_dy[2] = level_dy[1]*2

    # allocate matrix for storage at fine level
    # residual at fine level is already defined at global level
    prol_fine = zeros(Float64, level_nx[1]+1, level_ny[1]+1)

    # allocate matrix for storage at coarse levels
    source_coarse = zeros(Float64, level_nx[2]+1, level_ny[2]+1)
    u_numerical_coarse = zeros(Float64, level_nx[2]+1, level_ny[2]+1)

    # start main iteration loop
    for iteration_count = 1:maximum_iterations
        # call relaxation on fine grid and compute the numerical solution
        # for fixed number of iterations
        # @enter relax_multigrid(level_nx[1], level_ny[1], dx, dy, source, u_numerical, lambda,
        #                 relaxation_f2c)

        relax_multigrid(level_nx[1], level_ny[1], dx, dy, source, u_numerical, lambda,
                        tiny, relaxation_f2c)

        # check for convergence only for finest grid
        # compute the residual and L2 norm

        compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, residual)
        # write results only for finest residual
        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
                break
        end

        # restrict the residual from fine level to coarse level
        # @enter restriction(level_nx[1], level_ny[1], level_nx[2], level_ny[2], residual, source_coarse)

        restriction(level_nx[1], level_ny[1], level_nx[2], level_ny[2], residual, source_coarse)

        # set solution zero on coarse grid
        u_numerical_coarse[:,:] = zeros(level_nx[2]+1, level_ny[2]+1)

        # solve on the coarsest level and relax V3 times
        relax_multigrid(level_nx[2], level_ny[2], level_dx[2], level_dy[2],
                        source_coarse, u_numerical_coarse, lambda, tiny, relaxation_coarsest)

        # prolongate solution from coarse level to fine level
        prolongation(level_nx[2], level_ny[2], level_nx[1], level_ny[1], u_numerical_coarse, prol_fine)

        # correct the solution on fine level
        for j = 2:level_nx[1] for i = 2:level_ny[1]
                u_numerical[i,j] = u_numerical[i,j] + prol_fine[i,j]
        end end

        # relax v2 times
        relax_multigrid(level_nx[1], level_ny[1], dx, dy, source, u_numerical, lambda,
                        tiny, relaxation_c2f)
    end
end
