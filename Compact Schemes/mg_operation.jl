#------------------------- Restriction -----------------------------------------
# Arguments:
# nx_fine, ny_fine = number of grid points on fine level
# nx_coarse, ny_coarse = number of grid points on coarse level
# residual = residual at fine level
# source_coarse = source term at coarse level (∇2ϕ = ϵ_f)
#-------------------------------------------------------------------------------
function restriction(nx_fine, ny_fine, nx_coarse, ny_coarse, residual, source_coarse)

    for j = 2:ny_coarse for i = 2:nx_coarse
        # grid index for fine grid for the same coarse point
        center = 4.0*residual[2*i-1, 2*j-1]
        # E, W, N, S with respect to coarse grid point in fine grid
        grid = 2.0*(residual[2*i-1, 2*j-1+1] + residual[2*i-1, 2*j-1-1] +
                    residual[2*i-1+1, 2*j-1] + residual[2*i-1-1, 2*j-1])
        # NE, NW, SE, SW with respect to coarse grid point in fine grid
        corner = 1.0*(residual[2*i-1+1, 2*j-1+1] + residual[2*i-1+1, 2*j-1-1] +
                      residual[2*i-1-1, 2*j-1+1] + residual[2*i-1-1, 2*j-1-1])
        # restriction using trapezoidal rule
        source_coarse[i,j] = (center + grid + corner)/16.0

    end end

    # restriction for boundary points bottom and top
    for j = 1:nx_coarse+1
        # bottom boundary i = 1
        source_coarse[1,j] = residual[1, 2*j-1]
        # top boundary i = ny_coarse+1
        source_coarse[ny_coarse+1] = residual[ny_fine+1, 2*j-1]
    end

    # restriction for boundary poinys left and right
    for i = 1:ny_coarse+1
        # left boundary j = 1
        source_coarse[i,1] = residual[2*i-1,1]
        # right boundary nx_coarse+1
        source_coarse[i,nx_coarse+1] = residual[2*i-1, ny_fine+1]
    end
end

#-------------------------- Prolongation ---------------------------------------
# Arguments:
# nx_coarse, ny_coarse = number of grid points on coarse level
# nx_fine, ny_fine = number of grid points on fine level
# u_numerical_coarse = numerical solution at coarse level
# prol_fine = correction at fine level
#-------------------------------------------------------------------------------
function prolongation(nx_coarse, ny_coarse, nx_fine, ny_fine, u_numerical_coarse, prol_fine)

    for j = 1:ny_coarse for i = 1:nx_coarse
        # direct injection at center point
        prol_fine[2*i-1, 2*j-1] = u_numerical_coarse[i,j]
        # east neighnour on fine grid corresponding to coarse grid point
        prol_fine[2*i-1, 2*j-1+1] = 0.5*(u_numerical_coarse[i,j] + u_numerical[i,j+1])
        # north neighbout on fine grid corresponding to coarse grid point
        prol_fine[2*i-1+1, 2*j-1] = 0.5*(u_numerical_coarse[i,j] + u_numerical_coarse[i+1,j])
        # NE neighbour on fine grid corresponding to coarse grid point
        prol_fine[2*i-1+1, 2*j-1+1] = 0.25*(u_numerical_coarse[i,j] + u_numerical_coarse[i,j+1] +
                                            u_numerical_coarse[i+1,j] + u_numerical_coarse[i+1,j+1])

    end end

    # update boundary points
    for i = 1:ny_coarse+1
        # left boundary j = 1
        prol_fine[2*i-1,1] = u_numerical_coarse[i,1]
        # right boundary j = nx_fine+1
        prol_fine[2*i-1, nx_fine+1] = u_numerical_coarse[i,nx_coarse+1]
    end

    for j = 1:nx_coarse+1
        #bottom boundary i = 1
        prol_fine[1,2*j-1] = u_numerical_coarse[1,j]
        # top boundary i =  ny_fine+1
        prol_fine[ny_fine+1,2*j-1] = u_numerical_coarse[ny_coarse+1,j]
    end
end
