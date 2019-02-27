include("problem_assignment.jl")
include("residualcalculation.jl")
include("gauss_seidel.jl")
include("solver.jl")

using CPUTime
using Printf

# read data from text file for input parameters
file_input = open("input.txt")
input_lines = readlines(file_input)
input_parameters = Array{Float64}(undef, 15)
counter = 1
for line in input_lines
#    m = match(r"^([0-9]+)", line)
    m = split(line, "\t")
    input_parameters[counter] = parse(Float64, m[1])
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
tolerance               = Float64(input_parameters[7])
omega                   = Float64(input_parameters[8])
relaxation_f2c          = Int32(input_parameters[9])
relaxation_c2f          = Int32(input_parameters[10])
relaxation_coarsest     = Int32(input_parameters[11])
maximum_iterations      = Int32(input_parameters[12])
tiny                    = Float64(input_parameters[13])
lambda                  = Float64(input_parameters[14])
flag_order              = Int32(input_parameters[15])
k = 1

# create output file for L2-norm
output = open("output.txt", "w");

# create text file for initial and final field
field_initial = open("field_initial.txt", "w");
field_final = open("field_final.txt", "w");

write(field_initial, "variables =\"x\",\"y\",\"f\",\"u\",\"ue\" \n")
write(field_initial, "zone f=point i = ", string(nx), ",j = ", string(ny), "\n")

write(field_final, "variables =\"x\",\"y\",\"f\",\"u\",\"ue\", \"e\" \n")
write(field_final, "zone f=point i = ", string(nx), ",j = ", string(ny), "\n")

# Assign the domain size based on initial problem
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

# calculation of grid spacing
dx = (x_right - x_left)/nx
dy = (y_top - y_bottom)/ny

# allocate array for x and y position of grids, exact solution and source term
x_position      = Array{Float64}(undef, nx+1)
y_position      = Array{Float64}(undef, ny+1)
u_exact         = Array{Float64}(undef, nx+1, ny+1)
source          = Array{Float64}(undef, nx+1, ny+1)
u_numerical     = Array{Float64}(undef, nx+1, ny+1)

# assign the computational domain, exact solution, source term ,
# initial condition, boundary condition
calculate_grid_position(x_position, y_position, nx, ny, dx, dy, x_left,
                        x_right, y_top, y_bottom)

# calculate the source term and exact solution based on the problem flag
assign_problem(nx, ny, x_position, y_position, source, u_exact, flag_problem,
               lambda)

# assign initial condition (zero or random) based on start flag
initial_condition(u_numerical, nx, ny, flag_start)

# assign boundary conditions (exact solution) for the domain
boundary_condition(u_numerical, u_exact, nx, ny)

# calculate residual and initial l2 norm of the residual
residual         = zeros(Float64, nx+1, ny+1)
initial_rms      = 0.0
rms              = 0.0
for j = 1:ny+1 for i = 1:nx+1
    write(field_initial, @sprintf("%.16f",x_position[i])," ", @sprintf("%.16f", y_position[j]), " ",
          @sprintf("%.16f", source[i,j])," ", @sprintf("%.16f", u_numerical[i,j])," ",
          @sprintf("%.16f", u_exact[i,j]), " \n")
end end

@time begin
# call the solver function to calculate the numerical solution
    solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
           initial_rms, maximum_iterations, tiny, lambda, output, flag_solver)

end
for j = 1:ny+1 for i = 1:nx+1
    write(field_final, @sprintf("%.16f",x_position[i])," ", @sprintf("%.16f", y_position[j]), " ",
          @sprintf("%.16f", source[i,j])," ", @sprintf("%.16f", u_numerical[i,j])," ",
          @sprintf("%.16f", u_exact[i,j])," ", @sprintf("%.16f",(u_numerical[i,j]-u_exact[i,j]))," \n")
end end

close(field_initial)
close(field_final)
close(output);
