include("residualcalculation.jl")
include("relax_multigrid.jl")
include("mg_operation.jl")

#-------------------------Multigrid Solver with 2 levels------------------------
# this function executes multigrid framework with two levels
#
#-------------------------------------------------------------------------------
function mg3(dx, dy, nx, ny, residual, source, u_numerical, rms,
            initial_rms, maximum_iterations, tiny, lambda, output, n_level)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    # define 3D matrix for u_multigrid
    u_multigrid = Matrix{Float64}[]
    source_multigrid = Matrix{Float64}[]
    # fine mesh numerical solution and source temp at first level of 3D matrix
    push!(u_multigrid, u_numerical)
    push!(source_multigrid, source)

    # compute initial residual
    compute_residual(nx, ny, dx, dy, source_multigrid[1], u_multigrid[1], residual, lambda)
    # compute initial L-2 norm
    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    println("0", " ", rms, " ", rms/initial_rms)

    if nx < (2^n_level)
        println("Number of levels exceeds the possible number.\n")
    end
    #allocate memory for grid size at different levels
    level_nx = zeros(Int64, n_level)
    level_ny = zeros(Int64, n_level)
    level_dx = zeros(Float64, n_level)
    level_dy = zeros(Float64, n_level)

    # initialize the mesh details at fine level
    level_nx[1] = nx
    level_ny[1] = ny
    level_dx[1] = dx
    level_dy[1] = dy
    # calculate mesh details for coarse levels and allocate matirx for
    # numerical solution and error restreicted from upper level
    println(level_nx)
    for i = 2:n_level
        level_nx[i] = Int64(level_nx[i-1]/2)
        level_ny[i] = Int64(level_ny[i-1]/2)
        level_dx[i] = level_dx[i-1]*2
        level_dy[i] = level_dy[i-1]*2
        println(i, " ", level_nx)
        # allocate matrix for storage at coarse levels
        source_coarse = zeros(Float64, level_nx[i]+1, level_ny[i]+1)
        u_numerical_coarse = zeros(Float64, level_nx[i]+1, level_ny[i]+1)

        push!(u_multigrid, u_numerical_coarse)
        push!(source_multigrid, source_coarse)
    end

    # allocate matrix for storage at fine level
    # residual at fine level is already defined at global level
    prol_fine = zeros(Float64, level_nx[1]+1, level_ny[1]+1)
    # temporaty residual which is restricted to coarse mesh error
    # the size keeps on changing

    # start main iteration loop
    for iteration_count = 1:maximum_iterations
        # call relaxation on fine grid and compute the numerical solution
        # for fixed number of iterations
        # @enter relax_multigrid(level_nx[1], level_ny[1], dx, dy, source, u_numerical, lambda,
        #                 relaxation_f2c)

        relax_multigrid(level_nx[1], level_ny[1], level_dx[1], level_dy[1],
                        source_multigrid[1], u_multigrid[1], lambda, tiny, relaxation_f2c)

        # println(u_multigrid[1])
        # check for convergence only for finest grid
        # compute the residual and L2 norm

        compute_residual(nx, ny, dx, dy, source_multigrid[1], u_multigrid[1], residual, lambda)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, residual)
        # write results only for finest residual
        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
                break
        end
        # from second level to coarsest level

        # for second level temporary residual is taken from fine mesh level
            temp_residual = zeros(Float64, level_nx[1]+1, level_ny[1]+1)
            temp_residual = residual

            # restrict reisudal from (k-1)th level to kth level
            restriction(level_nx[1], level_ny[1], level_nx[2], level_ny[2], temp_residual,
                        source_multigrid[2])

            # solution at kth level to zero
            u_multigrid[2][:,:] = zeros(level_nx[2]+1, level_ny[2]+1)

            relax_multigrid(level_nx[2], level_ny[2], level_dx[2], level_dy[2],
                        source_multigrid[2], u_multigrid[2], lambda, tiny, relaxation_f2c)
            # println(u_multigrid[2])

            temp_residual = zeros(Float64, level_nx[2]+1, level_ny[2]+1)
            compute_residual(level_nx[2], level_ny[2], level_dx[2], level_dy[2],
                             source_multigrid[2], u_multigrid[2], temp_residual, lambda)

            restriction(level_nx[2], level_ny[2], level_nx[3], level_ny[3], temp_residual,
                        source_multigrid[3])

            # solution at kth level to zero
            u_multigrid[3][:,:] = zeros(level_nx[3]+1, level_ny[3]+1)

            relax_multigrid(level_nx[3], level_ny[3], level_dx[3], level_dy[3],
                            source_multigrid[3], u_multigrid[3], lambda, tiny, relaxation_coarsest)

            # println(u_multigrid[3])

            prol_fine = zeros(Float64, level_nx[2]+1, level_ny[2]+1)
            prolongation(level_nx[3], level_ny[3], level_nx[2], level_ny[2],
                         u_multigrid[3], prol_fine)

            # println(prol_fine)
            for j = 2:level_nx[2] for i = 2:level_ny[2]
                    u_multigrid[2][i,j] = u_multigrid[2][i,j] + prol_fine[i,j]
            end end
            # println("check")
            # println(u_multigrid[2])
            # println(source_multigrid[2])

            relax_multigrid(level_nx[2], level_ny[2], level_dx[2], level_dy[2],
                            source_multigrid[2], u_multigrid[2], lambda, tiny, relaxation_c2f)

            # println(u_multigrid[2])
            prol_fine = zeros(Float64, level_nx[1]+1, level_ny[1]+1)
            prolongation(level_nx[2], level_ny[2], level_nx[1], level_ny[1],
                         u_multigrid[2], prol_fine)

            for j = 2:level_nx[1] for i = 2:level_ny[1]
                    u_multigrid[1][i,j] = u_multigrid[1][i,j] + prol_fine[i,j]
            end end

            relax_multigrid(level_nx[1], level_ny[1], level_dx[1], level_dy[1], source_multigrid[1],
                            u_multigrid[1], lambda, tiny, relaxation_c2f)
            # println(u_multigrid[1])
    end
    u_numerical = u_multigrid[1]
end
