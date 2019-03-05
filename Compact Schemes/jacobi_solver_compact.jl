include("residualcalculation.jl")
# include("residualcalculation_compact.jl")

#-------------------------------------------------------------------------------
# This function performs the jacobi iteration to compute the numerical
# solution at every step. Numerical solution is updated after residual is
# calculated over the whole domain
#
# Arguments:
#   nx, ny = number of grid in x and y direction
#   dx, dy = grid spacing in x and y direction
#   source = source term on rhs of poisson equation ∇2ϕ = S
#   u_numerical = numerical solution for different problems
#   residual = residual at each grid location
#   rms, initial_rms = root mean squire of residual
#   maximum_iterations = maximum number of iterations
#
# Algorithm:
# 10    r^(k+1) = S - ∇^2(ϕ^k)
# 20    ϕ^(k+1) = ϕ^k + ωr^(k+1)
# 30    calculate residual rms for ϕ^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function jacobi_solver_compact(dx, dy, nx, ny, residual, source, u_numerical, rms,
    initial_rms, maximum_iterations, lambda, output)
    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm_compact(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", rms, " ", rms/initial_rms)

    factor = -12.0/(5*dx*dx) - 12.0/(5*dy*dy) - lambda*lambda

    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = -12.0/(50.0*dx*dx) + 6.0/(5.0*dy*dy)
    ne = nw = se = sw = 6/(50*dx*dx) + 6/(50*dy*dy)
    cc = 12/(5*dx*dx) + 12/(5*dy*dy)
    lambda2 = lambda*lambda
    # coefficients for source term
    beta    = 1/10.0
    beta2   = 1/100.0

    for iteration_count = 1:maximum_iterations

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for i = 2:nx for j = 2:ny
        # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            f_grid = (source[i+1,j] + source[i-1,j] + source[i,j+1] + source[i,j-1])*beta
        # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            f_corner = (source[i+1,j+1] + source[i+1,j-1] + source[i-1,j+1] + source[i-1,j-1])*beta2

        #total source term for compact scheme
            F = source[i,j] + f_grid + f_corner

        # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*u_numerical[i+1,j] + ww*u_numerical[i-1,j] +
                     nn*u_numerical[i,j+1] + ss*u_numerical[i,j-1]
        # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*u_numerical[i+1,j+1] + nw*u_numerical[i-1,j+1] +
                       se*u_numerical[i+1,j-1] + sw*u_numerical[i-1,j-1]

            X = x_grid + x_corner

        # calculate residual
            residual[i,j] = F + lambda2*u_numerical[i,j] + cc*u_numerical[i,j] - X

        end end

        for i = 2:nx for j = 2:ny
            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor
        end end

        compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

        # compute the l2norm of residual
        rms = compute_l2norm_compact(nx, ny, residual)

        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
            break
        end
    end
    max_error = maximum(abs.(residual))
    write(output, "L-2 Norm = ", string(rms), " \n");
    write(output, "Maximum Norm = ", string(max_error), " \n");
    write(output, "Iterations = ", string(count), " \n");
    close(residual_plot)
end

#-------------------------------------------------------------------------------
# Relaxation using Jacobi method with 4th order compact scheme
# in multigrid framework
#-------------------------------------------------------------------------------
function jacobi_solver_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)

    residual = zeros(Float64, nx+1, ny+1)
    factor = -12.0/(5.0*dx*dx) - 12.0/(5.0*dy*dy) - lambda*lambda
    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = -12.0/(50.0*dx*dx) + 6.0/(5.0*dy*dy)
    ne = nw = se = sw = 6.0/(50.0*dx*dx) + 6.0/(50.0*dy*dy)
    cc = 12.0/(5.0*dx*dx) + 12.0/(5.0*dy*dy)
    lambda2 = lambda*lambda

    # coefficients for source term
    beta    = 1/10.0
    beta2   = 1/100.0
    for iteration_count = 1:V

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for j = 2:ny for i = 2:nx
        # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            f_grid = (source[i+1,j] + source[i-1,j] + source[i,j+1] + source[i,j-1])*beta
        # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            f_corner = (source[i+1,j+1] + source[i+1,j-1] + source[i-1,j+1] + source[i-1,j-1])*beta2

        #total source term for compact scheme
            F = source[i,j] + f_grid + f_corner

        # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*u_numerical[i+1,j] + ww*u_numerical[i-1,j] +
                     nn*u_numerical[i,j+1] + ss*u_numerical[i,j-1]
        # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*u_numerical[i+1,j+1] + nw*u_numerical[i-1,j+1] +
                       se*u_numerical[i+1,j-1] + sw*u_numerical[i-1,j-1]

            X = x_grid + x_corner

        # calculate residual
            residual[i,j] = F + lambda2*u_numerical[i,j] + cc*u_numerical[i,j] - X
        end end

        for j = 2:ny for i = 2:nx
        # update numerical solution
            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor
        end end
    end
end
