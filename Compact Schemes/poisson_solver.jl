include("./problem_assignment.jl")
# using problem_assignment

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

# close the file
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

x_left = 0.0
x_right = 0.0
y_top = 0.0
y_bottom = 0.0
# Assign the domain size based on initial problem
#assign_domain(x_left, x_right, y_top, y_bottom, flag_problem)

#function assign_domain(x_left::Float64, x_right::Float64, y_top::Float64, y_bottom::Float64, flag_problem::Int32)
    if flag_problem == 1
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
#end

# calculation of grid spacing
dx = (x_right - x_left)/nx
dy = (y_top - y_bottom)/ny

# allocate array for x and y position of grids, exact solution and source term
x_position      = Array{Float32}(undef, nx+1)
y_position      = Array{Float32}(undef, ny+1)
u_exact         = Array{Float32}(undef, nx+1, ny+1)
source          = Array{Float32}(undef, nx+1, ny+1)
u_numerical     = Array{Float32}(undef, nx+1, ny+1)

# assign the computational domain
calculate_grid_position(x_position, y_position, nx, ny, dx, dy, x_left,
                        x_right, y_top, y_bottom)

assign_problem(x_position, y_position, source, u_exact, flag_problem)

initial_condition(u_numerical, nx, ny, flag_start)
