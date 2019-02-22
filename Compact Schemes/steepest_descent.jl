include("residualcalculation.jl")

function steepest_descent(dx, dy, nx, ny, residual, source, u_numerical, rms,
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
# 00    r^0 = S - ∇^2(ϕ^0)
# 10    d = ∇^2(r^k)
# 20    aa = <r,r>
# 30    bb = <d,r>
# 40    cc = aa/bb
# 50    ϕ^(k+1) = ϕ^k + cc*r^k
# 60    r^(k+1) = r^k - cc*r^k
# 30    calculate rms for r^(k+1) and go to 10 if rms < tolerance
#-------------------------------------------------------------------------------
    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms
    # println(initial_rms)

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

        println(iteration_count, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
            break
        end
    end
end
