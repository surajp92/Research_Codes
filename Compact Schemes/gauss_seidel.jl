include("residualcalculation.jl")

function gauss_seidel(dx, dy, nx, ny, residual, source, u_numerical, rms,
                      initial_rms, maximum_iterations)
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
    compute_residual(nx, ny, dx, dy, source, u_numerical, residual)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms

    factor = -2.0/dx^2 - 2.0/dy^2
    for iteration_count = 1:maximum_iterations

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for j = 2:ny for i = 2:nx
            residual[i,j] = source[i,j] - (u_numerical[i+1,j] - 2*u_numerical[i,j] +
                            u_numerical[i-1,j])/dx^2 - (u_numerical[i,j+1] -
                            2*u_numerical[i,j] + u_numerical[i,j-1])/dy^2

            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor
        end end

        compute_residual(nx, ny, dx, dy, source, u_numerical, residual)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, residual)

        println(iteration_count, " ", rms/initial_rms)

        if (rms/initial_rms) <= tolerance
            break
        end
    end
end
