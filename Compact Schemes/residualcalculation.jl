
function compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)
#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#    dx, dy = grid spacing in x and y direction
#    u_numerical = numericla solution
#
# Returns/ assigns:
#   residual = residual between source term and numericla solution at each grid
#-------------------------------------------------------------------------------

# residual = f + λ^2u - ∇^2u
    for i = 2:nx for j = 2:ny
        residual[i,j] = source[i,j] + lambda*lambda*u_numerical[i,j] -
                    (u_numerical[i+1,j] - 2*u_numerical[i,j] + u_numerical[i-1,j])/dx^2 -
                    (u_numerical[i,j+1] - 2*u_numerical[i,j] + u_numerical[i,j-1])/dy^2
    end end
end

function compute_l2norm(nx, ny, residual)
#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#   residual = residual between source term and numericla solution at each grid
#
# Retrns/ assigns:
#   rms =  root mean square of residual
#-------------------------------------------------------------------------------
    rms = 0.0
    for i = 2:nx for j = 2:ny
        rms = rms + residual[i,j]^2
    end end

    rms = sqrt(rms/(nx-1)*(ny-1))
    return rms
end
