include("problem_assignment.jl")
include("residualcalculation.jl")
include("gauss_seidel.jl")
include("solver.jl")
#import .pa

using Profile
using ProfileView
using BenchmarkTools

# clearconsole()

using CPUTime
using Printf
using ASTInterpreter2

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

flag = zeros(Int64,6)
relaxcount = zeros(Int64, 3)

# asssign input parameters
nx                      = Int32(input_parameters[1])
ny                      = Int32(input_parameters[1])
flag[1]                 = Int32(input_parameters[2]) # flag for the problem
flag[2]                 = Int32(input_parameters[3]) # flag for the solver
flag[3]                 = Int32(input_parameters[4]) # flag for multigrid
flag[4]                 = Int32(input_parameters[5]) # flag for staring condition
flag[5]                 = Int32(input_parameters[6]) # flag for restriction operator
tolerance               = Float64(input_parameters[7])
omega                   = Float64(input_parameters[8])
relaxcount[1]           = Int32(input_parameters[9]) # number of iterations from fine to coarse
relaxcount[2]           = Int32(input_parameters[10]) # number of iterations from coarse to fine
relaxcount[3]           = Int32(input_parameters[11]) # number of iterations on coarsest
maximum_iterations      = Int32(input_parameters[12])
tiny                    = Float64(input_parameters[13])
lambda                  = Float64(input_parameters[14])
flag[6]                 = Int32(input_parameters[15]) # flag for order of accuracy
k = 1


# create output file for L2-norm
output = open("output.txt", "w");
write(output, "Residual details: \n");
# create text file for initial and final field
field_initial = open("field_initial.txt", "w");
field_final = open("field_final.txt", "w");

write(field_initial, "variables =\"x\",\"y\",\"f\",\"u\",\"ue\" \n")
write(field_initial, "zone f=point i = ", string(nx), ",j = ", string(ny), "\n")

write(field_final, "variables =\"x\",\"y\",\"f\",\"u\",\"ue\", \"e\" \n")
write(field_final, "zone f=point i = ", string(nx), ",j = ", string(ny), "\n")


# Assign the domain size based on initial problem
    if flag[1] == 1
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
# function calculate_grid_position(x_position, y_position, nx, ny, dx, dy, x_left,
#                         x_right, y_top, y_bottom)
#     # assign x position for each grid point
#     for i = 1:nx+1
#         x_position[i] = x_left + dx*(i-1)
#     end
#
#     # assign x position for each grid point
#     for i = 1:ny+1
#         y_position[i] = y_bottom + dy*(i-1)
#     end
# end

calculate_grid_position(x_position, y_position, nx, ny, dx, dy, x_left,
                        x_right, y_top, y_bottom)


# function assign_problem(nx, ny, x_position, y_position, source, u_exact,
#                         flag_problem, lambda)
#         if flag_problem == 1
#         # a test case from Moin's book "Engineering NUmerical Analysis"
#             for i = 1:ny+1 for j = 1:nx+1
#
#                 source[i,j]  = -2.0 * (2.0 - x_position[i]^2 - y_position[j]^2)
#
#                 u_exact[i,j] = (x_position[i]^2 - 1.0)*(y_position[j]^2 - 1.0)
#             end end
#
#         elseif flag_problem == 2
#         # a test case from Moin's book "Engineering NUmerical Analysis"
#             c1 = (1.0/16.0)^2
#             c2 = -2.0*pi*pi
#             for i = 1:nx+1 for j = 1:ny+1
#
#                 source[i,j] = c2 * sin(pi * x_position[i]) * sin(pi * y_position[j]) +
#                               c2 * sin(16.0 * pi * x_position[i]) * sin(16.0 * pi * y_position[j])
#
#                 u_exact[i,j] = sin(pi * x_position[i]) * sin(pi * y_position[j]) +
#                                c1 * sin(16.0 * pi * x_position[i]) * sin(16.0 * pi * y_position[j])
#             end end
#
#         elseif flag_problem == 3
#         # a test case for sin function
#             c2 = -2.0*pi*pi
#             for i = 1:nx+1 for j = 1:ny+1
#
#                 source[i,j] = c2 * sin(pi*x_position[i]) * sin(pi*y_position[j])
#
#                 u_exact[i,j] = sin(pi*x_position[i]) * sin(pi*y_position[j])
#
#             end end
#         elseif flag_problem == 4
#         # exponential function
#             for i = 1:nx+1 for j = 1:ny+1
#
#                 source[i,j] = (x_position[i]^2 + y_position[j]^2) *
#                               exp(x_position[i] * y_position[j])
#
#                 u_exact[i,j] = exp(x_position[i] * y_position[j])
#             end end
#         else
#         # Helmholtz function with k = 1, l = 1
#             k = 0.5
#             l = 0.5
#             c = -(lambda^2 + (k*k + l*l)*pi^2)
#             for i = 1:nx+1 for j = 1:ny+1
#
#                 source[i,j] = c*sin(pi*k*x_position[i])*sin(pi*l*y_position[j])
#
#                 u_exact[i,j] = sin(pi*k*x_position[i]) * sin(pi*l*y_position[j])
#             end end
#         end
# end

# calculate the source term and exact solution based on the problem flag
assign_problem(nx, ny, x_position, y_position, source, u_exact, flag[1],
               lambda)

# function initial_condition(u_numerical, nx, ny, flag_start)
#     if flag_start == 1
#         for i = 1:nx+1 for j = 1:ny+1
#             u_numerical[i,j] = 0.0
#         end end
#     else
#         for i = 1:nx+1 for j = 1:ny+1
#             u_numerical[i,j] = rand()
#         end end
#     end
# end

# assign initial condition (zero or random) based on start flag
initial_condition(u_numerical, nx, ny, flag[4])

# function boundary_condition(u_numerical, u_exact, nx, ny)
#     for i = 1:nx+1
#         u_numerical[i,1] = u_exact[i,1]
#         u_numerical[i, ny+1] = u_exact[i, ny+1]
#     end
#
#     for j = 1:ny+1
#         u_numerical[1,j] = u_exact[1,j]
#         u_numerical[nx+1,j] = u_exact[nx+1,j]
#     end
# end

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

# @enter solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
#       initial_rms, maximum_iterations, tiny, lambda, output, flag_solver)

# @profile solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
#        initial_rms, maximum_iterations, tiny, lambda, output, flag, relaxcount)

val, t, bytes, gctime, memallocs = @timed begin
# @time begin
# call the solver function to calculate the numerical solution
    # println("a")

    solver(dx, dy, nx, ny, residual, source, u_numerical, rms,
           initial_rms, maximum_iterations, tiny, lambda, output, flag, relaxcount,
           omega)

end

# ProfileView.view()
# Profile.print()

u_error = zeros(nx+1, ny+1)
rms_error = 0.0

for j = 1:ny+1 for i = 1:nx+1
    u_error[i,j] = u_numerical[i,j] - u_exact[i,j]
end end

rms_error = compute_l2norm(nx, ny, u_error)
max_error = maximum(abs.(u_error))

println("Error details:");
println("L-2 Norm = ", rms_error);
println("Maximum Norm = ", max_error);
print("CPU Time = ", t);

write(output, "Error details: \n");
write(output, "L-2 Norm = ", string(rms_error), " \n");
write(output, "Maximum Norm = ", string(max_error), " \n");
write(output, "CPU Time = ", string(t), " \n");

for j = 1:ny+1 for i = 1:nx+1
    write(field_final, @sprintf("%.16f",x_position[i])," ", @sprintf("%.16f", y_position[j]), " ",
          @sprintf("%.16f", source[i,j])," ", @sprintf("%.16f", u_numerical[i,j])," ",
          @sprintf("%.16f", u_exact[i,j])," ", @sprintf("%.16f",(u_numerical[i,j]-u_exact[i,j]))," \n")
end end

close(field_initial)
close(field_final)
close(output);
