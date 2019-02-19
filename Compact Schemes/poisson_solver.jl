
# read data from text file for input parameters
file_input = open("input.txt")
input_lines = readlines(file_input)
input_parameters = Array{Float32}(undef, 13)
counter = 1
for line in input_lines
    m = match(r"^([0-9]+)", line)
    input_parameters[counter] = parse(Float32, m[1])
    global counter += 1
end

close(file_input)

# asssign input parameters
nx                      = Int32(input_parameters[1])
ny                      = Int32(input_parameters[1])
flag_problem            = Int32(input_parameters[2])
flag_solver             = Int32(input_parameters[3])
flag_multigrid          = Int32(input_parameters[4])
flag_start              = Int32(input_parameters[5])
flag_restriction        = Int32(input_parameters[6])
tolerance               = Float32(input_parameters[7])
omega                   = Float32(input_parameters[8])
relaxation_f2c          = Int32(input_parameters[9])
relaxaation_c2f         = Int32(input_parameters[10])
relazation_coarsest     = Int32(input_parameters[11])
maximum_iterations      = Int32(input_parameters[12])
tiny                    = Float32(input_parameters[13])
k = 1

# Assign the domain size based on initial problem
if problem == 1
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
