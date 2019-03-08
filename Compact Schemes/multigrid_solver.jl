include("residualcalculation.jl")
# include("residualcalculation_compact.jl")
include("relax_multigrid.jl")
include("mg_operation.jl")

#-------------------------Multigrid Solver with 2 levels------------------------
# this function executes multigrid framework with two levels
#
#-------------------------------------------------------------------------------
function multigrid_solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
            initial_rms, maximum_iterations, tiny, lambda, output, n_level, relaxcount, flag,
            omega)

    relaxation_f2c          = relaxcount[1]
    relaxation_c2f          = relaxcount[2]
    relaxation_coarsest     = relaxcount[3]
    flag_solver             = flag[2]
    flag_order              = flag[6]
    # println(flag_solver)
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
    function_residual(nx, ny, dx, dy, source_multigrid[1], u_multigrid[1], residual, lambda, flag_order)

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

    for i = 2:n_level
        level_nx[i] = Int64(level_nx[i-1]/2)
        level_ny[i] = Int64(level_ny[i-1]/2)
        level_dx[i] = level_dx[i-1]*2
        level_dy[i] = level_dy[i-1]*2

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
    temp_residual = zeros(Float64, level_nx[1]+1, level_ny[1]+1)
    # start main iteration loop
    for iteration_count = 1:maximum_iterations
        # call relaxation on fine grid and compute the numerical solution
        # for fixed number of iterations
        # @enter relax_multigrid(level_nx[1], level_ny[1], dx, dy, source, u_numerical, lambda,
        #                 relaxation_f2c)

        relax_multigrid(level_nx[1], level_ny[1], level_dx[1], level_dy[1],
                        source_multigrid[1], u_multigrid[1], lambda, tiny,
                        relaxation_f2c, flag_solver, flag_order, omega)
        # println(u_multigrid[1])
        # check for convergence only for finest grid
        # compute the residual and L2 norm
        function_residual(level_nx[1], level_ny[1], level_dx[1], level_dy[1],
                source_multigrid[1], u_multigrid[1], residual, lambda, flag_order)

        # compute the l2norm of residual
        rms = compute_l2norm(level_nx[1], level_ny[1], residual)
        # write results only for finest residual
        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
                break
        end
        # from second level to coarsest level
        for k = 2:n_level
            if k == 2
                # for second level temporary residual is taken from fine mesh level
                temp_residual = residual
            else
                # from third level onwards residual is computed for (k-1) level
                # which will be restricted to kth lvel error
                temp_residual = zeros(Float64, level_nx[k-1]+1, level_ny[k-1]+1)
                function_residual(level_nx[k-1], level_ny[k-1], level_dx[k-1], level_dy[k-1],
                        source_multigrid[k-1], u_multigrid[k-1], temp_residual, lambda, flag_order)
            end
            # restrict reisudal from (k-1)th level to kth level
            restriction(level_nx[k-1], level_ny[k-1], level_nx[k], level_ny[k], temp_residual,
                        source_multigrid[k])

            # solution at kth level to zero
            u_multigrid[k][:,:] = zeros(level_nx[k]+1, level_ny[k]+1)

            # solve (∇^-λ^2)ϕ = ϵ on coarse grid (kthe level)
            if k < n_level
                relax_multigrid(level_nx[k], level_ny[k], level_dx[k], level_dy[k],
                            source_multigrid[k], u_multigrid[k], lambda, tiny,
                            relaxation_f2c, flag_solver, flag_order, omega)
            elseif k == n_level
                relax_multigrid(level_nx[k], level_ny[k], level_dx[k], level_dy[k],
                            source_multigrid[k], u_multigrid[k], lambda, tiny,
                            relaxation_coarsest, flag_solver, flag_order, omega)
            end
            # println(u_multigrid[k])
        end

        for k = n_level:-1:2
            # temporary matrix for correction storage at (k-1) th level
            # solution prolongated from kth level to (k-1)th level
            prol_fine = zeros(Float64, level_nx[k-1]+1, level_ny[k-1]+1)

            # prolongate solution from (k)th level to (k-1)th level
            prolongation(level_nx[k], level_ny[k], level_nx[k-1], level_ny[k-1],
                         u_multigrid[k], prol_fine)

            for j = 2:level_nx[k-1] for i = 2:level_ny[k-1]
                    u_multigrid[k-1][i,j] = u_multigrid[k-1][i,j] + prol_fine[i,j]
            end end

            relax_multigrid(level_nx[k-1], level_ny[k-1], level_dx[k-1], level_dy[k-1],
                            source_multigrid[k-1], u_multigrid[k-1], lambda, tiny,
                            relaxation_c2f, flag_solver, flag_order, omega)

            #println(u_multigrid[k-1])
        end
    end
    u_numerical = u_multigrid[1]

    max_residual = maximum(abs.(residual))
    write(output, "L-2 Norm = ", string(rms), " \n");
    write(output, "Maximum Norm = ", string(max_residual), " \n");
    write(output, "Iterations = ", string(count), " \n");
    close(residual_plot)
end
