include("residualcalculation.jl")

#------------------------------ Steepest Descent -------------------------------
# This function performs the steepest descent iteration to compute the numerical
# solution at every step. Numerical solution is updated while the residuals
# are being calculated
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
# 00    r^0 = S - ∇^2(ϕ^0)
# 10    d = ∇^2(r^k)
# 20    aa = <r,r>
# 30    bb = <d,r>
# 40    cc = aa/bb
# 50    ϕ^(k+1) = ϕ^k + cc*r^k
# 60    r^(k+1) = r^k - cc*r^k
# 30    calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function steepest_descent(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", rms, " ", rms/initial_rms)

    del_residual    = zeros(Float64, nx+1, ny+1)

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx

            del_residual[i,j] = (residual[i+1,j] - 2*residual[i,j] + residual[i-1,j])/(dx^2) +
                                (residual[i,j+1] - 2*residual[i,j] + residual[i,j-1])/(dy^2) -
                                lambda*lambda*residual[i,j]
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_residual[i,j]*residual[i,j]
        end end
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of residual
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * residual[i,j]
        end end

        # update the residual by removiing some component of previous residual
        for j = 2:ny for i = 2:nx
            residual[i,j] = residual[i,j] - cc * del_residual[i,j]
        end end

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, residual)
        # println(rms)

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

#------------------------------- Steepest Descent multigrid --------------------
# This function performs steepst descent algorithm for fixed number of iterations
# during relaxation of solution in multigrid framework
#-------------------------------------------------------------------------------
function steepest_descent_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)

    # allocate temporary residual matrix
    residual = zeros(Float64, nx+1, ny+1)
    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    del_residual    = zeros(Float64, nx+1, ny+1)

    # start calculation
    for iteration_count = 1:V

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx
            del_residual[i,j] = (residual[i+1,j] - 2*residual[i,j] + residual[i-1,j])/(dx^2) +
                                (residual[i,j+1] - 2*residual[i,j] + residual[i,j-1])/(dy^2) -
                                lambda*lambda*residual[i,j]
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_residual[i,j]*residual[i,j]
        end end
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of residual
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * residual[i,j]
        end end

        # update the residual by removiing some component of previous residual
        for j = 2:ny for i = 2:nx
            residual[i,j] = residual[i,j] - cc * del_residual[i,j]
        end end
    end
end

#************************** Compact Scheme  ************************************

#-------------------------------------------------------------------------------
# This function performs the gauss seidel iteration to compute the numerical
# solution at every step. Numerical solution is updated while the residuals
# are being calculated
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
# 00    r^0 = S - ∇^2(ϕ^0)
# 10    d = ∇^2(r^k)
# 20    aa = <r,r>
# 30    bb = <d,r>
# 40    cc = aa/bb
# 50    ϕ^(k+1) = ϕ^k + cc*r^k
# 60    r^(k+1) = r^k - cc*r^k
# 30    calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function steepest_descent_compact(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm_compact(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", rms, " ", rms/initial_rms)

    del_residual    = zeros(Float64, nx+1, ny+1)

    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = 6.0/(5.0*dy*dy) - 12.0/(50.0*dx*dx)
    ne = nw = se = sw = 6.0/(50.0*dx*dx) + 6.0/(50.0*dy*dy)
    cc1 = 12.0/(5.0*dx*dx) + 12.0/(5.0*dy*dy)
    lambda2 = lambda*lambda

    # source term scaling factor
    gg = 100.0/144.0

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate (∇^2-λ^2)(residual)
        for j = 2:ny for i = 2:nx
            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*residual[i+1,j] + ww*residual[i-1,j] +
                     nn*residual[i,j+1] + ss*residual[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*residual[i+1,j+1] + nw*residual[i-1,j+1] +
                       se*residual[i+1,j-1] + sw*residual[i-1,j-1]
            X = x_grid + x_corner

            del_residual[i,j] = (X - cc1*residual[i,j] -lambda2*residual[i,j])*gg
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_residual[i,j]*residual[i,j]
        end end
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of residual
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * residual[i,j]
        end end

        # update the residual by removiing some component of previous residual
        for j = 2:ny for i = 2:nx
            residual[i,j] = residual[i,j] - cc * del_residual[i,j]
        end end

        # compute the l2norm of residual
        rms = compute_l2norm_compact(nx, ny, residual)
        # println(rms)

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

#------------------------------- Steepest Descent multigrid --------------------
# This function performs steepst descent algorithm with 4th order
# compact scheme for fixed number of iterations during relaxation of solution in
# multigrid framework
#-------------------------------------------------------------------------------
function steepest_descent_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)

    residual = zeros(Float64, nx+1, ny+1)
    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    del_residual    = zeros(Float64, nx+1, ny+1)

    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = 6.0/(5.0*dy*dy) - 12.0/(50.0*dx*dx)
    ne = nw = se = sw = 6.0/(50.0*dx*dx) + 6.0/(50.0*dy*dy)
    cc1 = 12.0/(5.0*dx*dx) + 12.0/(5.0*dy*dy)
    lambda2 = lambda*lambda

    # source term scaling factor
    gg = 100.0/144.0

    # start calculation
    for iteration_count = 1:V

        # calculate (∇^2-λ^2)(residual)
        for j = 2:ny for i = 2:nx
            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*residual[i+1,j] + ww*residual[i-1,j] +
                     nn*residual[i,j+1] + ss*residual[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*residual[i+1,j+1] + nw*residual[i-1,j+1] +
                       se*residual[i+1,j-1] + sw*residual[i-1,j-1]
            X = x_grid + x_corner

            del_residual[i,j] = (X - cc1*residual[i,j] -lambda2*residual[i,j])*gg
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_residual[i,j]*residual[i,j]
        end end
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of residual
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * residual[i,j]
        end end

        # update the residual by removiing some component of previous residual
        for j = 2:ny for i = 2:nx
            residual[i,j] = residual[i,j] - cc * del_residual[i,j]
        end end
    end
end
