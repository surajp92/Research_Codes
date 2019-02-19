module problem_assignment


function assign_problem(x_position, y_position, nx, ny)
    # assign x position for each grid point
    for i = 1:nx+1
        x_position[i] = x_left + dx*(i-1)
    end

    # assign x position for each grid point
    for i = 1:ny+1
        y_position[i] = y_bottom + dy*(i-1)
    end
end

end  # module problem_assignment
