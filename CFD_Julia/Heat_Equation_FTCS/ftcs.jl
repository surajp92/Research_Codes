clearconsole()

using CPUTime
using Printf
using Plots

x_l = -1.0
x_r = 1.0
dx = 0.025
nx = Int64((x_r-x_l)/dx)

dt = 0.0025
t = 1.0
nt = Int64(t/dt)

α = 1/(pi*pi)

x = Array{Float64}(undef, nx+1)
u_e = Array{Float64}(undef, nx+1)
u_n = Array{Float64}(undef, nt+1, nx+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)  # location of each grid point
    u_n[1,i] = -sin(pi*x[i]) # initial condition @ t=0
    u_e[i] = -exp(-t)*sin(pi*x[i]) # initial condition @ t=0
end

u_n[1,1] = 0.0
u_n[1,nx+1] = 0.0

α1 = α*dt/(dx*dx)

for k = 2:nt+1
    for i = 2:nx
        u_n[k,i] = α1*(u_n[k-1,i+1] - 2.0*u_n[k-1,i] + u_n[k-1,i-1]) +
                   u_n[k-1,i]
    end
    u_n[k,1] = 0.0
    u_n[k,nx+1] = 0.0
end

p1 = plot(x,u_e,lw = 4,xlabel="X", color = :red, ylabel = "U", xlims=(minimum(x),maximum(x)),
     grid=(:none), label = "Exact")

p2 = plot(x,u_n[nt+1,:],lw = 4,xlabel="X", color = :blue, ylabel = "U", xlims=(minimum(x),maximum(x)),
     grid=(:none), label = "t=1")

plot(p1, p2)
