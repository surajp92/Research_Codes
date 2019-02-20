#module problem_assignment

function calculate_grid_position(x_position, y_position, nx, ny, dx, dy, x_left,
                                x_right, y_top, y_bottom)
#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#    dx, dy = grid spacing in x and y direction
#    x_left, x_right, y_top, y_bottom = boindary locations
#
# Returns/ assigns
#    x_position = an array with x-coordinates of each grid location
#    y_position = an array with y-coordinates of each grid location
#-------------------------------------------------------------------------------
    # assign x position for each grid point
    for i = 1:nx+1
        x_position[i] = x_left + dx*(i-1)
    end

    # assign x position for each grid point
    for i = 1:ny+1
        y_position[i] = y_bottom + dy*(i-1)
    end
end

function assign_problem(x_position, y_position, source, u_exact, flag_problem)
#-------------------------------------------------------------------------------
# Arguments:
#    x_position = an array with x-coordinates of each grid location
#    y_position = an array with y-coordinates of each grid location
#    flag_problem = flag for the problem [1]- Moin, [2]- harmonic(16), [3]- sin
#                   [4]- exponential
#
# Returns/ assigns
#    source = source term on rhs of poisson equation ∇2ϕ = S
#    u_exact = exact solution for different problems
#-------------------------------------------------------------------------------
    if flag_problem == 1
    # a test case from Moin's book "Engineering NUmerical Analysis"
        for i = 1:ny+1 for j = 1:nx+1

            source[i,j]  = -2.0 * (2.0 - x_position[i]^2 - y_position[j]^2)

            u_exact[i,j] = (x_position[i]^2 - 1.0)*(y_position[j]^2 - 1.0)
        end end

    elseif flag_problem == 2
    # a test case from Moin's book "Engineering NUmerical Analysis"
        c1 = (1.0/16.0)^2
        c2 = -2.0*pi*2
        for i = 1:nx+1 for j = 1:ny+1

            source[i,j] = c2 * sin(pi * x_position[i]) * sin(pi * y_position[j])
                        + c2 * sin(16.0 * pi * x_position[i]) *
                        sin(16.0 * pi * y_position[j])

            u_exact[i,j] = sin(pi * x_position[i]) * sin(pi * y_position[j]) +
                           c1 * sin(16.0 * pi * x_position[i]) *
                           sin(16.0 * pi * y_position[j])
        end end

    elseif flag_problem == 3
    # a test case for sin function
        c2 = -2.0*pi*pi
        for i = 1:nx+1 for j = 1:ny+1

            source[i,j] = c2 * sin(pi * x_position[i]) * sin(pi * y_position[j])

            u_exact[i,j] = sin(pi * x_position[i]) * sin(pi * y_position[j])

        end end
    else
    # exponential function
        for i = 1:nx+1 for j = 1:ny+1

            source[i,j] = (x_position[i]^2 + y_position[j]^2) *
                          exp(x_position[i] * y_position[j])

            u_exact[i,j] = exp(x_position[i] * y_position[j])
        end end
    end
end
#end  # module problem_assignment

function initial_condition(u_numerical, nx, ny, flag_start)
#-------------------------------------------------------------------------------
# Arguments:
#    nx, ny = number of grid in x and y direction
#    flag_start = flag for initial condition [1]- zero, [2]- random
#
# Returns/ assigns:
#   u_numerical = numerical solution
#-------------------------------------------------------------------------------
    if flag_start == 1
        for i = 1:nx+1 for j = 1:ny+1
            u_numerical[i,j] = 0.0
        end end
    else
        for i = 1:nx+1 for j = 1:ny+1
            u_numerical[i,j] = rand()
        end end
    end
end

function boundary_condition(u_numerical, u_exact, nx, ny)
#-------------------------------------------------------------------------------
# Arguments:
#    nx,ny = number of grid in x and y direction
#    u_numerical = numerical solution
#    u_exact = exact solution
#
# Returns/ assigns:
#    boundary points of u_numerical
#-------------------------------------------------------------------------------

    for i = 1:nx+1
        u_numerical[i,1] = u_exact[i,1]
        u_numerical[i, ny+1] = u_exact[i, ny+1]
    end

    for j = 1:ny+1
        u_numerical[1,j] = u_exact[1,j]
        u_numerical[nx+1,j] = u_exact[nx+1,j]
    end
end
