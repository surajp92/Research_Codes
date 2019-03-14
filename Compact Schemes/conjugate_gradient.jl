include("residualcalculation.jl")

#--------------------------- Conjugate gradient algorithm-----------------------
# This function performs the conjugate gradient algorithm to compute the numerical
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
# 00    r^0 = S - ∇^2(ϕ^0); p^0 = r^0
# 10    d^(k+1) = ∇^2(p^k)
# 20    aa = <r,r>
# 30    bb = <d,p>
# 40    cc = aa/bb
# 50    ϕ^(k+1) = ϕ^k + cc*p^k
# 60    bb = <r,r> or (bb = aa)
# 70    r^(k+1) = r^k -cc*d^k
# 60    aa = <r^(k+1),r^(k+1) >
# 70    cc = aa/bb
# 80    p^(k+1) = r^(k+1) - cc*p^k
# 30    calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function conjugate_gradient(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.plt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    write(residual_plot, "zone T=\"", string(nx), " x ", string(ny), "\"\n")

    count = 0.0

    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    println(initial_rms)
    # allocate the matric for direction and set the initial direction (conjugate vector)
    p = zeros(Float64, nx+1, ny+1)

    # asssign conjugate vector to initial residual
    for j = 1:ny+1 for i = 1:nx+1
        p[i,j] = residual[i,j]
    end end

    del_p    = zeros(Float64, nx+1, ny+1)

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx
            del_p[i,j] = (p[i+1,j] - 2*p[i,j] + p[i-1,j])/(dx^2) +
                         (p[i,j+1] - 2*p[i,j] + p[i,j-1])/(dy^2) -
                         lambda*lambda*p[i,j]
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_p[i,j]*p[i,j]
        end end
        # cc = <r,r>/<d,p>
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of conjugate vector
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * p[i,j]
        end end

        # bb = <r,r> = aa (calculated in previous loop)
        bb = aa
        aa = 0.0

        # update the residual by removing some component of previous residual
        for j = 1:ny for i = 1:nx
            residual[i,j] = residual[i,j] - cc * del_p[i,j]
            aa = aa + residual[i,j]*residual[i,j]
        end end
        # cc = <r-cd, r-cd>/<r,r>
        cc = aa/bb

        # update the conjugate vector
        for j = 1:ny for i = 1:nx
            p[i,j] = residual[i,j] + cc * p[i,j]
        end end

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

#------------------------------- Conjugate gradient multigrid ------------------
# This function performs the conjugate gradient algorithm  during relaxation
# for each level for Fixed numbe rof times
#-------------------------------------------------------------------------------
function conjugate_gradient_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)

    # allocate temporary residual matrix
    residual = zeros(Float64, nx+1, ny+1)
    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    p = zeros(Float64, nx+1, ny+1)

    # asssign conjugate vector to initial residual
    for j = 1:ny+1 for i = 1:nx+1
        p[i,j] = residual[i,j]
    end end

    del_p    = zeros(Float64, nx+1, ny+1)

    # start calculation
    for iteration_count = 1:V

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx
            del_p[i,j] = (p[i+1,j] - 2*p[i,j] + p[i-1,j])/(dx^2) +
                         (p[i,j+1] - 2*p[i,j] + p[i,j-1])/(dy^2) -
                         lambda*lambda*p[i,j]
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_p[i,j]*p[i,j]
        end end
        # cc = <r,r>/<d,p>
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of conjugate vector
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * p[i,j]
        end end

        # bb = <r,r> = aa (calculated in previous loop)
        bb = aa
        aa = 0.0

        # update the residual by removing some component of previous residual
        for j = 1:ny for i = 1:nx
            residual[i,j] = residual[i,j] - cc * del_p[i,j]
            aa = aa + residual[i,j]*residual[i,j]
        end end
        # cc = <r-cd, r-cd>/<r,r>
        cc = aa/(bb+tiny)

        # update the conjugate vector
        for j = 1:ny for i = 1:nx
            p[i,j] = residual[i,j] + cc * p[i,j]
        end end
    end
    # println("Relaxation")
    # println(u_numerical)
end

#**************************** Compact scheme ***********************************
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
# 00    r^0 = S - ∇^2(ϕ^0); p^0 = r^0
# 10    d^(k+1) = ∇^2(p^k)
# 20    aa = <r,r>
# 30    bb = <d,p>
# 40    cc = aa/bb
# 50    ϕ^(k+1) = ϕ^k + cc*p^k
# 60    bb = <r,r> or (bb = aa)
# 70    r^(k+1) = r^k -cc*d^k
# 60    aa = <r^(k+1),r^(k+1) >
# 70    cc = aa/bb
# 80    p^(k+1) = r^(k+1) - cc*p^k
# 30    calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function conjugate_gradient_compact(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.plt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    write(residual_plot, "zone T=\"", string(nx), " x ", string(ny), "\"\n")

    count = 0.0
    println("CG Compact")
    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm_compact(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", initial_rms, " ", rms/initial_rms)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    p = zeros(Float64, nx+1, ny+1)

    # asssign conjugate vector to initial residual
    for j = 1:ny+1 for i = 1:nx+1
        p[i,j] = residual[i,j]
    end end

    del_p    = zeros(Float64, nx+1, ny+1)

    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = 6.0/(5.0*dy*dy) - 12.0/(50.0*dx*dx)
    ne = nw = se = sw = 6.0/(50.0*dx*dx) + 6.0/(50.0*dy*dy)
    cc1 = 12.0/(5.0*dx*dx) + 12.0/(5.0*dy*dy)
    lambda2 = lambda*lambda

    # source term scaling factor
    beta    = 1/10.0
    beta2   = 1/100.0
    gg = 100.0/144.0

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*p[i+1,j] + ww*p[i-1,j] +
                     nn*p[i,j+1] + ss*p[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*p[i+1,j+1] + nw*p[i-1,j+1] +
                       se*p[i+1,j-1] + sw*p[i-1,j-1]
            X = x_grid + x_corner

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            p_grid = (p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1])*beta
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            p_corner = (p[i+1,j+1] + p[i+1,j-1] + p[i-1,j+1] + p[i-1,j-1])*beta2

            #total source term for compact scheme
            P = p[i,j] + p_grid + p_corner

            # del_p[i,j] = (X - cc1*p[i,j] - lambda2*p[i,j])*gg
            del_p[i,j] = (X - cc1*p[i,j] - lambda2*P)*gg

        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_p[i,j]*p[i,j]
        end end
        # cc = <r,r>/<d,p>
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of conjugate vector
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * p[i,j]
        end end

        # bb = <r,r> = aa (calculated in previous loop)
        bb = aa
        aa = 0.0
        # bb = 0.0

        # update the residual by removing some component of previous residual
        for j = 2:ny for i = 2:nx
            # bb = bb + residual[i,j]*residual[i,j]
            residual[i,j] = residual[i,j] - cc * del_p[i,j]
            aa = aa + residual[i,j]*residual[i,j]
        end end
        # cc = <r-cd, r-cd>/<r,r>
        cc = aa/(bb+tiny)

        # update the conjugate vector
        for j = 2:ny for i = 2:nx
            p[i,j] = residual[i,j] + cc * p[i,j]
        end end

        # compute the l2norm of residual
        rms = compute_l2norm_compact(nx, ny, residual)

        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
            break
        end
    end
    max_residual = maximum(abs.(residual))
    write(output, "L-2 Norm = ", string(rms), " \n");
    write(output, "Maximum Norm = ", string(max_residual), " \n");
    write(output, "Iterations = ", string(count), " \n");
    close(residual_plot)
end


#------------------------------- Conjugate gradient multigrid ------------------
# This function performs the conjugate gradient algorithm with 4th order
# compact schemer  during relaxation for each level for Fixed numbe rof times
#-------------------------------------------------------------------------------
function conjugate_gradient_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)

    # allocate temporary residual matrix
    residual = zeros(Float64, nx+1, ny+1)
    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    p = zeros(Float64, nx+1, ny+1)

    # asssign conjugate vector to initial residual
    for j = 1:ny+1 for i = 1:nx+1
        p[i,j] = residual[i,j]
    end end

    del_p    = zeros(Float64, nx+1, ny+1)

    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = 6.0/(5.0*dy*dy) - 12.0/(50.0*dx*dx)
    ne = nw = se = sw = 6.0/(50.0*dx*dx) + 6.0/(50.0*dy*dy)
    cc1 = 12.0/(5.0*dx*dx) + 12.0/(5.0*dy*dy)
    lambda2 = lambda*lambda

    # source term scaling factor
    beta    = 1/10.0
    beta2   = 1/100.0
    gg = 100.0/144.0

    # start calculation
    for iteration_count = 1:V

        # calculate ∇^2(residual)
        for j = 2:ny for i = 2:nx
            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*p[i+1,j] + ww*p[i-1,j] +
                     nn*p[i,j+1] + ss*p[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*p[i+1,j+1] + nw*p[i-1,j+1] +
                       se*p[i+1,j-1] + sw*p[i-1,j-1]
            X = x_grid + x_corner

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            p_grid = (p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1])*beta
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            p_corner = (p[i+1,j+1] + p[i+1,j-1] + p[i-1,j+1] + p[i-1,j-1])*beta2

            #total source term for compact scheme
            P = p[i,j] + p_grid + p_corner

            # del_p[i,j] = (X - cc1*p[i,j] - lambda2*p[i,j])*gg
            del_p[i,j] = (X - cc1*p[i,j] - lambda2*P)*gg
        end end

        aa = 0.0
        bb = 0.0
        # calculate aa, bb, cc. cc is the distance parameter(α_n)
        for j = 2:ny for i = 2:nx
            aa = aa + residual[i,j]*residual[i,j]
            bb = bb + del_p[i,j]*p[i,j]
        end end
        # cc = <r,r>/<d,p>
        cc = aa/(bb + tiny)

        # update the numerical solution by adding some component of conjugate vector
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + cc * p[i,j]
        end end

        # bb = <r,r> = aa (calculated in previous loop)
        bb = aa
        aa = 0.0

        # update the residual by removing some component of previous residual
        for j = 1:ny for i = 1:nx
            residual[i,j] = residual[i,j] - cc * del_p[i,j]
            aa = aa + residual[i,j]*residual[i,j]
        end end
        # cc = <r-cd, r-cd>/<r,r>
        cc = aa/bb

        # update the conjugate vector
        for j = 1:ny for i = 1:nx
            p[i,j] = residual[i,j] + cc * p[i,j]
        end end
    end
end
