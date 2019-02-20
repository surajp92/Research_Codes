include("residualcalculation.jl")

function jacobi_solver()
    factor = -2.0/dx^2 - 2.0/dy^2
    for iteration_count = 1:maximum_iterations

        for i = 2:nx for j = 2:ny
            residual[i,j] = source[i,j] - (u_numerical[i+1,j] - 2*u_numerical[i,j] +
                            u_numerical[i-1,j])/dx^2 - (u_numerical[i,j+1] -
                            2*u_numerical[i,j] + u_numerical[i,j-1])/dy^2

        end end

        compute_residual(nx, ny, dx, dy, source, u_numerical, residual)

        rms = compute_l2norm(nx, ny, residual)

        println(rms/initial_rms)

        if rms/initial_rms >= tolerance
            break
        end
    end
end
