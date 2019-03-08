#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
function function_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda, flag_order)
    if flag_order == 1
        compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)
    elseif flag_order == 2
        compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)
    end
end

#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#    dx, dy = grid spacing in x and y direction
#    u_numerical = numericla solution
#
# Returns/ assigns:
#   residual = residual between source term and numericla solution at each grid
# residual = f + λ^2u - ∇^2u
#-------------------------------------------------------------------------------
function compute_residual(nx, ny, dx, dy, source, u_numerical, residual, lambda)

    for j = 2:ny for i = 2:nx
        d2udx2 = (u_numerical[i+1,j] - 2*u_numerical[i,j] + u_numerical[i-1,j])/(dx^2)
        d2udy2 = (u_numerical[i,j+1] - 2*u_numerical[i,j] + u_numerical[i,j-1])/(dy^2)
        residual[i,j] = source[i,j]  - d2udx2 - d2udy2 + lambda*lambda*u_numerical[i,j]
    end end
end

#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#   residual = residual between source term and numericla solution at each grid
#
# Retrns/ assigns:
#   rms =  root mean square of residual
#-------------------------------------------------------------------------------
function compute_l2norm(nx, ny, residual)

    rms = 0.0
    # println(residual)
    for j = 2:ny for i = 2:nx
        rms = rms + residual[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx-1)*(ny-1)))
    return rms
end

#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#    dx, dy = grid spacing in x and y direction
#    u_numerical = numericla solution
#
# Returns/ assigns:
#   residual = residual between source term and numericla solution at each grid
# residual = f + λ^2u - ∇^2u
#-------------------------------------------------------------------------------
function compute_residual_compact(nx, ny, dx, dy, source, u_numerical, residual, lambda)


    # calculate constant coefficients
    ee = ww = 6.0/(5.0*dx*dx) - 12.0/(50.0*dy*dy)
    nn = ss = 6.0/(5.0*dy*dy) - 12.0/(50.0*dx*dx)
    ne = nw = se = sw = 6.0/(50.0*dx*dx) + 6.0/(50.0*dy*dy)
    cc = 12.0/(5.0*dx*dx) + 12.0/(5.0*dy*dy)
    lambda2 = lambda*lambda
    # coefficients for source term
    beta    = 1/10.0
    beta2   = 1/100.0

    gg = 100.0/144.0

    for i = 2:nx for j = 2:ny

        # alternate 4th order accurate method
        # F = 0.5*dx*dx*(source[i+1,j] + source[i-1,j] + source[i,j+1] + source[i,j-1] + 8*source[i,j])
        # x_grid = 4.0*(u_numerical[i+1,j] + u_numerical[i-1,j] + u_numerical[i,j+1] + u_numerical[i,j-1])
        # x_corner = (u_numerical[i+1,j+1] + u_numerical[i-1,j+1] + u_numerical[i+1,j-1] + u_numerical[i-1,j-1])
        # residual[i,j] = (source[i,j] - x_grid - x_corner + 20*u_numerical[i,j])

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
        residual[i,j] = (F + lambda2*u_numerical[i,j] + cc*u_numerical[i,j] - X)*gg
    end end
end

#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#   residual = residual between source term and numericla solution at each grid
#
# Retrns/ assigns:
#   rms =  root mean square of residual
#-------------------------------------------------------------------------------
function compute_l2norm_compact(nx, ny, residual)
    rms = 0.0
    # println(residual)
    for i = 2:nx for j = 2:ny
        rms = rms + residual[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx-1)*(ny-1)))
    return rms
end
