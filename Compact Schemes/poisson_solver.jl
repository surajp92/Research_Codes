# using DelimitedFiles
# file_input = open("input.txt")
# nx = parse(Int,split(readline(file_input), "#"))
# close(file_input)
nx = 128
ny = 128
initial_problem = 1

# Assign the domain size based on initial problem
if initial_problem == 1
    x_left      = -1.0
    x_right     = 1.0

    y_bottom    = -1.0
    y_top       = 1.0
else
    x_left      = 0.0
    x_right     = 1.0

    y_bottom    = 0.0
    y_top       = 1.0
end

# calculation of grid spacing
dx = (x_right - x_left)/nx
dy = (y_top - y_bottom)/ny

# calculation of spatial coordinates
x_position = Array{Float32}(undef, nx+1)
y_position = Array{Float32}(undef, ny+1)

# assign x position for each grid point
for i = 1:nx+1
    x_position[i] = x_left + dx*(i-1)
end

# assign x position for each grid point
for i = 1:ny+1
    y_position[i] = y_bottom + dy*(i-1)
end

#
