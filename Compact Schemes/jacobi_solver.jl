include("residualcalculation.jl")

#-------------------------Jacobi Solver-----------------------------------------
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
function jacobi_solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
                       initial_rms, maximum_iterations, lambda, output)

    # create text file for writing residual history
    residual_plot = open("residual.txt", "w")
    write(residual_plot, "variables =\"k\",\"rms\",\"rms/rms0\"\n")
    count = 0.0

    compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    rms = compute_l2norm(nx, ny, residual)

    initial_rms = rms

    factor = -2.0/dx^2 - 2.0/dy^2 - lambda*lambda
    for iteration_count = 1:maximum_iterations

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for j = 2:nx for i = 2:ny
            residual[i,j] = source[i,j] + lambda*lambda*u_numerical[i,j]-
                        (u_numerical[i+1,j] - 2*u_numerical[i,j] + u_numerical[i-1,j])/dx^2 -
                        (u_numerical[i,j+1] - 2*u_numerical[i,j] + u_numerical[i,j-1])/dy^2

        end end
        for j = 2:nx for i = 2:ny
            u_numerical[i,j] = u_numerical[i,j] + omega * residual[i,j]/factor
        end end

        compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

        # compute the l2norm of residual
        rms = compute_l2norm(nx, ny, residual)

        write(residual_plot, string(iteration_count), " ",string(rms), " ", string(rms/initial_rms)," \n");
        count = iteration_count

        println(iteration_count, " ", rms/initial_rms)

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


#-------------------------Jacobi Solver Multigrid-------------------------------
function jacobi_solver_mg(nx, ny, dx, dy, source, u_numerical, lambda, V)

    temp_residual = zeros(Float64, nx+1, ny+1)
    factor = -2.0/dx^2 - 2.0/dy^2 - lambda*lambda

    for iteration_count = 1:V

        # compute solution at next time step ϕ^(k+1) = ϕ^k + ωr^(k+1)
        for j = 2:nx for i = 2:ny
            temp_residual[i,j] = source[i,j] + lambda*lambda*u_numerical[i,j]-
                        (u_numerical[i+1,j] - 2*u_numerical[i,j] + u_numerical[i-1,j])/dx^2 -
                        (u_numerical[i,j+1] - 2*u_numerical[i,j] + u_numerical[i,j-1])/dy^2

        end end

        for j = 2:nx for i = 2:ny
            u_numerical[i,j] = u_numerical[i,j] + omega * temp_residual[i,j]/factor
        end end
    end
end
