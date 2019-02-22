include("residualcalculation.jl")

function conjugate_gradient(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations, tiny, lambda)
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
                         lambda*lambda*residual[i,j]
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
        for j = 2:ny-1 for i = 2:nx-1
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

        println(iteration_count, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
            break
        end
    end
end
