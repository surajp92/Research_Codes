include("residualcalculation.jl")

#-------------------------------------------------------------------------------
# This function performs the gauss seidel iteration to compute the numerical
# solution at every step. Numerical solution is updated while the residuals
# are being calculated
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
function gauss_seidel(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations, lambda, output, omega)
    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", initial_rms, " ", rms/initial_rms)

    factor = -2.0/dx^2 - 2.0/dy^2 - lambda*lambda
    for iteration_count = 1:maximum_iterations

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        # residual = f + λ^2u - ∇^2u
        for j = 2:ny for i = 2:nx
            d2udx2 = (u_numerical[i+1,j] - 2*u_numerical[i,j] + u_numerical[i-1,j])/(dx^2)
            d2udy2 = (u_numerical[i,j+1] - 2*u_numerical[i,j] + u_numerical[i,j-1])/(dy^2)
            residual[i,j] = source[i,j] + lambda*lambda*u_numerical[i,j] - d2udx2 - d2udy2

            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor
        end end

        compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, residual)

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


#-------------------------Gauss Seidel Multigrid-------------------------------
# This function performs the iteration of gauss seidel solver for fixed numer V
# Arguments:
#   nx, ny = number of grid in x and y direction
#   dx, dy = grid spacing in x and y direction
#   source = source term on rhs of poisson equation ∇2ϕ = S
#   (error restricted from fine levels for coarse levels)
#   u_numerical = numerical solution for different problems
#
#-------------------------------------------------------------------------------
function gauss_seidel_mg(nx, ny, dx, dy, source, u_numerical, lambda, V, omega)

    residual = zeros(Float64, nx+1, ny+1)
    factor = -2.0/dx^2 - 2.0/dy^2 - lambda*lambda

    for iteration_count = 1:V
        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for j = 2:nx for i = 2:ny
            residual[i,j] = source[i,j] + lambda*lambda*u_numerical[i,j]-
                        (u_numerical[i+1,j] - 2*u_numerical[i,j] + u_numerical[i-1,j])/dx^2 -
                        (u_numerical[i,j+1] - 2*u_numerical[i,j] + u_numerical[i,j-1])/dy^2

            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor

        end end

    end
    # println("Relaxation")
    # println(u_numerical)
end

#**************************** Compact Scheme ***********************************

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
function gauss_seidel_compact(dx, dy, nx, ny, residual, source, u_numerical, rms,
                       initial_rms, maximum_iterations, lambda, output, omega)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm_compact(nx, ny, residual)

    initial_rms = rms

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

    for iteration_count = 1:maximum_iterations

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

        # caclulate numerical solution ϕ^(k+1) = ϕ^k + ωr^(k+1)
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
# Relaxation using Gauss- seidel method with 4th order compact scheme
# in multigrid framework
#-------------------------------------------------------------------------------
function gauss_seidel_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, V, omega)

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

        # update numerical solution
            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor
        end end
    end
end
