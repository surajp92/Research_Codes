include("residualcalculation.jl")

#------------------------ Bi-conjugate gradient Stabilized ---------------------
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
# 00    r^0 = S - ∇^2(ϕ^0); f^0 = r^0
# 10    p^0 = 0.0, q^0 = 0.0...............................vector
# 20    α^0 = 1.0, ω^0 = 1.0, ρ^0 = 1.0....................scalar
# 30    ρ^(k+1) = <f^0, r^k>...............................scalar
# 40    β^k = (α^k/β^k)(ρ^(k+1)/ρ^k).......................scalar
# 50    p^(k+1) = r^k + β^k(p^k - ω^k*q^k).................vector
# 60    q^(k+1) = ∇^2(p^(k+1)).............................vector
# 70    α^(k+1) = ρ^(k+1)/<f^0, q^(k+1)>...................scalar
# 80    s^k = r^k - α^(k+1)*q^(k+1)........................vector
# 90    t^k = ∇^2(s^k).....................................vector
# 100   ω^(k+1) = <s^k, t^k>/<t^k, t^k>....................scalar
# 110   u^(k+1) = u^k + α^(k+1)*p^(k+1) + ω^(k+1)*s^k......vector
# 120   r^(k+1) = s^k - ω^(k+1)*t^k........................vector
# 130   calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function biconjugate_gradient_stab(dx, dy, nx, ny, residual, source, u_numerical,
                                   rms, initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.plt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    write(residual_plot, "zone T=\"", string(nx), " x ", string(ny), "\"\n")
    count = 0.0

    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", initial_rms, " ", rms/initial_rms)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    f_initial = zeros(Float64, nx+1, ny+1)

    # asssign frozen initial residual
    for j = 1:ny+1 for i = 1:nx+1
        f_initial[i,j] = residual[i,j]
    end end

    # initialize bi-conjugate vectors
    p    = zeros(Float64, nx+1, ny+1)
    q    = zeros(Float64, nx+1, ny+1)
    s    = zeros(Float64, nx+1, ny+1)
    t    = zeros(Float64, nx+1, ny+1)

    # initialize constant scalars
    alfa  = 1.0
    omega = 1.0
    rho   = 1.0

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate ρ^(k+1)
        rho_new = 0.0

        for j = 2:ny for i = 2:nx
            rho_new = rho_new + f_initial[i,j]*residual[i,j]
        end end

        # calculate β^k, rho^k = rho^(k+1)
        beta = (alfa * rho_new)/(omega * rho + tiny)
        rho = rho_new

        # calculate p^(k+1) = r^k + β^k(p^k - ω^k*q^k)
        for j = 2:ny for i = 2:nx
            p[i,j] = residual[i,j] + beta * (p[i,j] - omega * q[i,j])
        end end

        # apply filtering operator q^(k+1) = ∇^2(p^(k+1))
        for j = 2:ny for i = 2:nx
            q[i,j] = (p[i+1,j] - 2*p[i,j] + p[i-1,j])/(dx^2) +
                     (p[i,j+1] - 2*p[i,j] + p[i,j-1])/(dy^2) -
                     lambda*lambda*p[i,j]
        end end
        temp = 0.0

        # calculate α^(k+1) = ρ^(k+1)/<f^0, q^(k+1)>
        for j = 2:ny for i = 2:nx
            temp = temp + f_initial[i,j]*q[i,j]
        end end

        alfa = rho/(temp + tiny)

        # calcualte s^k = r^k - α^(k+1)*q^(k+1)
        for j = 2:ny for i = 2:nx
            s[i,j] = residual[i,j] - alfa * q[i,j]
        end end

        # apply filtering operator t^k = ∇^2(s^k)
        for j = 2:ny for i = 2:nx
            t[i,j] = (s[i+1,j] - 2*s[i,j] + s[i-1,j])/(dx^2) +
                     (s[i,j+1] - 2*s[i,j] + s[i,j-1])/(dy^2) -
                     lambda*lambda*s[i,j]
        end end

        # clculate ω^(k+1) = <s^k, t^k>/<t^k, t^k>
        aa = 0.0
        bb = 0.0
        for j = 2:ny for i = 2:nx
            aa = aa + s[i,j]*t[i,j]
            bb = bb + t[i,j]*t[i,j]
        end end

        omega = aa/(bb + tiny)

        # calculate numerical solution u^(k+1) = u^k + α^(k+1)*p^(k+1) + ω^(k+1)*s^k
        # calculate new residual r^(k+1) = s^k - ω^(k+1)*t^k
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + alfa * p[i,j] + omega * s[i,j]
            residual[i,j] = s[i,j] - omega * t[i,j]
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

#------------------------------- Bi-conjugate gradient multigrid ---------------
# This function performs the biconjugate gradient algorithm  during relaxation
# for each level for Fixed numbe rof times
#-------------------------------------------------------------------------------

function biconjugate_gradient_stab_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)

    # allocate temporary residual matrix
    residual = zeros(Float64, nx+1, ny+1)
    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    f_initial = zeros(Float64, nx+1, ny+1)

    # asssign frozen initial residual
    for j = 1:ny+1 for i = 1:nx+1
        f_initial[i,j] = residual[i,j]
    end end

    # initialize bi-conjugate vectors
    p    = zeros(Float64, nx+1, ny+1)
    q    = zeros(Float64, nx+1, ny+1)
    s    = zeros(Float64, nx+1, ny+1)
    t    = zeros(Float64, nx+1, ny+1)

    # initialize constant scalars
    alfa  = 1.0
    omega = 1.0
    rho   = 1.0

    # start calculation
    for iteration_count = 1:V

        # calculate ρ^(k+1)
        rho_new = 0.0

        for j = 2:ny for i = 2:nx
            rho_new = rho_new + f_initial[i,j]*residual[i,j]
        end end

        # calculate β^k, rho^k = rho^(k+1)
        beta = (alfa * rho_new)/(omega * rho + tiny)
        rho = rho_new

        # calculate p^(k+1) = r^k + β^k(p^k - ω^k*q^k)
        for j = 2:ny for i = 2:nx
            p[i,j] = residual[i,j] + beta * (p[i,j] - omega * q[i,j])
        end end

        # apply filtering operator q^(k+1) = ∇^2(p^(k+1))
        for j = 2:ny for i = 2:nx
            q[i,j] = (p[i+1,j] - 2*p[i,j] + p[i-1,j])/(dx^2) +
                     (p[i,j+1] - 2*p[i,j] + p[i,j-1])/(dy^2) -
                     lambda*lambda*p[i,j]
        end end
        temp = 0.0

        # calculate α^(k+1) = ρ^(k+1)/<f^0, q^(k+1)>
        for j = 2:ny for i = 2:nx
            temp = temp + f_initial[i,j]*q[i,j]
        end end

        alfa = rho/(temp + tiny)

        # calcualte s^k = r^k - α^(k+1)*q^(k+1)
        for j = 2:ny for i = 2:nx
            s[i,j] = residual[i,j] - alfa * q[i,j]
        end end

        # apply filtering operator t^k = ∇^2(s^k)
        for j = 2:ny for i = 2:nx
            t[i,j] = (s[i+1,j] - 2*s[i,j] + s[i-1,j])/(dx^2) +
                     (s[i,j+1] - 2*s[i,j] + s[i,j-1])/(dy^2) -
                     lambda*lambda*s[i,j]
        end end

        # clculate ω^(k+1) = <s^k, t^k>/<t^k, t^k>
        aa = 0.0
        bb = 0.0
        for j = 2:ny for i = 2:nx
            aa = aa + s[i,j]*t[i,j]
            bb = bb + t[i,j]*t[i,j]
        end end

        omega = aa/(bb + tiny)

        # calculate numerical solution u^(k+1) = u^k + α^(k+1)*p^(k+1) + ω^(k+1)*s^k
        # calculate new residual r^(k+1) = s^k - ω^(k+1)*t^k
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + alfa * p[i,j] + omega * s[i,j]
            residual[i,j] = s[i,j] - omega * t[i,j]
        end end
    end
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
# 00    r^0 = S - ∇^2(ϕ^0); f^0 = r^0
# 10    p^0 = 0.0, q^0 = 0.0...............................vector
# 20    α^0 = 1.0, ω^0 = 1.0, ρ^0 = 1.0....................scalar
# 30    ρ^(k+1) = <f^0, r^k>...............................scalar
# 40    β^k = (α^k/β^k)(ρ^(k+1)/ρ^k).......................scalar
# 50    p^(k+1) = r^k + β^k(p^k - ω^k*q^k).................vector
# 60    q^(k+1) = ∇^2(p^(k+1)).............................vector
# 70    α^(k+1) = ρ^(k+1)/<f^0, q^(k+1)>...................scalar
# 80    s^k = r^k - α^(k+1)*q^(k+1)........................vector
# 90    t^k = ∇^2(s^k).....................................vector
# 100   ω^(k+1) = <s^k, t^k>/<t^k, t^k>....................scalar
# 110   u^(k+1) = u^k + α^(k+1)*p^(k+1) + ω^(k+1)*s^k......vector
# 120   r^(k+1) = s^k - ω^(k+1)*t^k........................vector
# 130   calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
function biconjugate_gradient_stab_compact(dx, dy, nx, ny, residual, source, u_numerical,
                                   rms, initial_rms, maximum_iterations, tiny, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.plt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    write(residual_plot, "zone T=\"", string(nx), " x ", string(ny), "\"\n")

    count = 0.0

    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm_compact(nx, ny, residual)

    initial_rms = rms
    iteration_count = 0
    println(iteration_count, " ", initial_rms, " ", rms/initial_rms)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    f_initial = zeros(Float64, nx+1, ny+1)

    # asssign frozen initial residual
    for j = 1:ny+1 for i = 1:nx+1
        f_initial[i,j] = residual[i,j]
    end end

    # initialize bi-conjugate vectors
    p    = zeros(Float64, nx+1, ny+1)
    q    = zeros(Float64, nx+1, ny+1)
    s    = zeros(Float64, nx+1, ny+1)
    t    = zeros(Float64, nx+1, ny+1)

    # calculate constant coefficients
    ee = ww = 6/(5*dx*dx) - 12/(50*dy*dy)
    nn = ss = 6/(5*dy*dy) - 12/(50*dx*dx)
    ne = nw = se = sw = 6/(50*dx*dx) + 6/(50*dy*dy)
    cc1 = 12/(5*dx*dx) + 12/(5*dy*dy)
    lambda2 = lambda*lambda

    # source term scaling factor
    beta_s    = 1/10.0
    beta2_s   = 1/100.0
    gg = 100.0/144.0

    # initialize constant scalars
    alfa  = 1.0
    omega = 1.0
    rho   = 1.0

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate ρ^(k+1)
        rho_new = 0.0

        for j = 2:ny for i = 2:nx
            rho_new = rho_new + f_initial[i,j]*residual[i,j]
        end end

        # calculate β^k, rho^k = rho^(k+1)
        beta = (alfa * rho_new)/(omega * rho + tiny)
        rho = rho_new

        # calculate p^(k+1) = r^k + β^k(p^k - ω^k*q^k)
        for j = 2:ny for i = 2:nx
            p[i,j] = residual[i,j] + beta * (p[i,j] - omega * q[i,j])
        end end

        # apply filtering operator q^(k+1) = ∇^2(p^(k+1))
        for j = 2:ny for i = 2:nx

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*p[i+1,j] + ww*p[i-1,j] +
                     nn*p[i,j+1] + ss*p[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*p[i+1,j+1] + nw*p[i-1,j+1] +
                       se*p[i+1,j-1] + sw*p[i-1,j-1]
            X = x_grid + x_corner

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            p_grid = (p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1])*beta_s
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            p_corner = (p[i+1,j+1] + p[i+1,j-1] + p[i-1,j+1] + p[i-1,j-1])*beta2_s

            #total source term for compact scheme
            P = p[i,j] + p_grid + p_corner

            # q[i,j] = (X - cc1*p[i,j] - lambda2*p[i,j])*gg
             q[i,j] = (X - cc1*p[i,j] - lambda2*P)*gg

        end end
        temp = 0.0

        # calculate α^(k+1) = ρ^(k+1)/<f^0, q^(k+1)>
        for j = 2:ny for i = 2:nx
            temp = temp + f_initial[i,j]*q[i,j]
        end end

        alfa = rho/(temp + tiny)

        # calcualte s^k = r^k - α^(k+1)*q^(k+1)
        for j = 2:ny for i = 2:nx
            s[i,j] = residual[i,j] - alfa * q[i,j]
        end end

        # apply filtering operator t^k = ∇^2(s^k)
        for j = 2:ny for i = 2:nx

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*s[i+1,j] + ww*s[i-1,j] +
                     nn*s[i,j+1] + ss*s[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*s[i+1,j+1] + nw*s[i-1,j+1] +
                       se*s[i+1,j-1] + sw*s[i-1,j-1]
            X = x_grid + x_corner

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            s_grid = (s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1])*beta_s
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            s_corner = (s[i+1,j+1] + s[i+1,j-1] + s[i-1,j+1] + s[i-1,j-1])*beta2_s

            #total source term for compact scheme
            S = s[i,j] + s_grid + s_corner

            #t[i,j] = (X - cc1*s[i,j] - lambda2*s[i,j])*gg
            t[i,j] = (X - cc1*s[i,j] - lambda2*S)*gg

        end end

        # clculate ω^(k+1) = <s^k, t^k>/<t^k, t^k>
        aa = 0.0
        bb = 0.0
        for j = 2:ny for i = 2:nx
            aa = aa + s[i,j]*t[i,j]
            bb = bb + t[i,j]*t[i,j]
        end end

        omega = aa/(bb + tiny)

        # calculate numerical solution u^(k+1) = u^k + α^(k+1)*p^(k+1) + ω^(k+1)*s^k
        # calculate new residual r^(k+1) = s^k - ω^(k+1)*t^k
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + alfa * p[i,j] + omega * s[i,j]
            residual[i,j] = s[i,j] - omega * t[i,j]
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
    max_error = maximum(abs.(residual))
    write(output, "L-2 Norm = ", string(rms), " \n");
    write(output, "Maximum Norm = ", string(max_error), " \n");
    write(output, "Iterations = ", string(count), " \n");
    close(residual_plot)
end

#------------------------------- Bi-conjugate gradient multigrid ---------------
# This function performs the biconjugate gradient algorithm  with 4th order
# compact scheme during relaxation for each level for Fixed numbe rof times
#-------------------------------------------------------------------------------
function biconjugate_gradient_stab_compact_mg(nx, ny, dx, dy, source, u_numerical, lambda, tiny, V)

    # allocate temporary residual matrix
    residual = zeros(Float64, nx+1, ny+1)
    compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    # allocate the matric for direction and set the initial direction (conjugate vector)
    f_initial = zeros(Float64, nx+1, ny+1)

    # asssign frozen initial residual
    for j = 1:ny+1 for i = 1:nx+1
        f_initial[i,j] = residual[i,j]
    end end

    # initialize bi-conjugate vectors
    p    = zeros(Float64, nx+1, ny+1)
    q    = zeros(Float64, nx+1, ny+1)
    s    = zeros(Float64, nx+1, ny+1)
    t    = zeros(Float64, nx+1, ny+1)

    # calculate constant coefficients
    ee = ww = 6/(5*dx*dx) - 12/(50*dy*dy)
    nn = ss = 6/(5*dy*dy) - 12/(50*dx*dx)
    ne = nw = se = sw = 6/(50*dx*dx) + 6/(50*dy*dy)
    cc1 = 12/(5*dx*dx) + 12/(5*dy*dy)
    lambda2 = lambda*lambda

    # source term scaling factor
    beta_s    = 1/10.0
    beta2_s   = 1/100.0
    gg = 100.0/144.0

    # initialize constant scalars
    alfa  = 1.0
    omega = 1.0
    rho   = 1.0

    # start calculation
    for iteration_count = 1:maximum_iterations

        # calculate ρ^(k+1)
        rho_new = 0.0

        for j = 2:ny for i = 2:nx
            rho_new = rho_new + f_initial[i,j]*residual[i,j]
        end end

        # calculate β^k, rho^k = rho^(k+1)
        beta = (alfa * rho_new)/(omega * rho + tiny)
        rho = rho_new

        # calculate p^(k+1) = r^k + β^k(p^k - ω^k*q^k)
        for j = 2:ny for i = 2:nx
            p[i,j] = residual[i,j] + beta * (p[i,j] - omega * q[i,j])
        end end

        # apply filtering operator q^(k+1) = ∇^2(p^(k+1))
        for j = 2:ny for i = 2:nx

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*p[i+1,j] + ww*p[i-1,j] +
                     nn*p[i,j+1] + ss*p[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*p[i+1,j+1] + nw*p[i-1,j+1] +
                       se*p[i+1,j-1] + sw*p[i-1,j-1]
            X = x_grid + x_corner

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            p_grid = (p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1])*beta_s
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            p_corner = (p[i+1,j+1] + p[i+1,j-1] + p[i-1,j+1] + p[i-1,j-1])*beta2_s

            #total source term for compact scheme
            P = p[i,j] + p_grid + p_corner

            # q[i,j] = (X - cc1*p[i,j] - lambda2*p[i,j])*gg
            q[i,j] = (X - cc1*p[i,j] - lambda2*P)*gg

        end end
        temp = 0.0

        # calculate α^(k+1) = ρ^(k+1)/<f^0, q^(k+1)>
        for j = 2:ny for i = 2:nx
            temp = temp + f_initial[i,j]*q[i,j]
        end end

        alfa = rho/(temp + tiny)

        # calcualte s^k = r^k - α^(k+1)*q^(k+1)
        for j = 2:ny for i = 2:nx
            s[i,j] = residual[i,j] - alfa * q[i,j]
        end end

        # apply filtering operator t^k = ∇^2(s^k)
        for j = 2:ny for i = 2:nx

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            x_grid = ee*s[i+1,j] + ww*s[i-1,j] +
                     nn*s[i,j+1] + ss*s[i,j-1]
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            x_corner = ne*s[i+1,j+1] + nw*s[i-1,j+1] +
                       se*s[i+1,j-1] + sw*s[i-1,j-1]
            X = x_grid + x_corner

            # stencil corresponding to (i+1,j) (i-1,j) (i,j+1) (i,j-1)
            s_grid = (s[i+1,j] + s[i-1,j] + s[i,j+1] + s[i,j-1])*beta_s
            # stencil corresponding to (i+1,j+1) (i+1,j-1) (i-1,j+1) (i-1,j-1)
            s_corner = (s[i+1,j+1] + s[i+1,j-1] + s[i-1,j+1] + s[i-1,j-1])*beta2_s

            #total source term for compact scheme
            S = s[i,j] + s_grid + s_corner

            # t[i,j] = (X - cc1*s[i,j] - lambda2*s[i,j])*gg
            t[i,j] = (X - cc1*s[i,j] - lambda2*S)*gg

        end end

        # clculate ω^(k+1) = <s^k, t^k>/<t^k, t^k>
        aa = 0.0
        bb = 0.0
        for j = 2:ny for i = 2:nx
            aa = aa + s[i,j]*t[i,j]
            bb = bb + t[i,j]*t[i,j]
        end end

        omega = aa/(bb + tiny)

        # calculate numerical solution u^(k+1) = u^k + α^(k+1)*p^(k+1) + ω^(k+1)*s^k
        # calculate new residual r^(k+1) = s^k - ω^(k+1)*t^k
        for j = 2:ny for i = 2:nx
            u_numerical[i,j] = u_numerical[i,j] + alfa * p[i,j] + omega * s[i,j]
            residual[i,j] = s[i,j] - omega * t[i,j]
        end end
    end
end
