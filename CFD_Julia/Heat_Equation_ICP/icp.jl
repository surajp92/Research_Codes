using CPUTime
using Printf
using Plots

function compute_l2norm(nx,r)

    rms = 0.0

    for i = 2:nx
        rms = rms + r[i]^2
    end

    rms = sqrt(rms/((nx-1)))
    return rms
end

function tdma(a,b,c,r,s,e)
    up = Array{Float64}(undef, nx+1)
    for i = s+1:e+1
        b[i] = b[i] - a[i]*(c[i-1]/b[i-1])
        r[i] = r[i] - a[i]*(r[i-1]/b[i-1])
    end

    up[e+1] = r[e+1]/b[e+1]

    for i = e:-1:s
        up[i] = (r[i] - c[i]*up[i+1])/b[i]
    end
    return up
end

x_l = -1.0
x_r = 1.0
dx = 0.025
nx = Int64((x_r-x_l)/dx)

dt = 0.0025
t = 1.0
nt = Int64(t/dt)

α = 1/(pi*pi)

x = Array{Float64}(undef, nx+1)
tlist = Array{Float64}(undef, nt+1)
u_e = Array{Float64}(undef, nx+1)
u_n = Array{Float64}(undef, nt+1, nx+1)
a = Array{Float64}(undef, nx+1)
b = Array{Float64}(undef, nx+1)
c = Array{Float64}(undef, nx+1)
r = Array{Float64}(undef, nx+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)  # location of each grid point
    u_n[1,i] = -sin(pi*x[i]) # initial condition @ t=0
    u_e[i] = -exp(-t)*sin(pi*x[i]) # initial condition @ t=0
end

u_n[1,1] = 0.0
u_n[1,nx+1] = 0.0

s = 1
e = nx

for k = 2:nt+1
    i = 1
    a[i] = 0.0
    b[i] = 1.0
    c[i] = 0.0
    r[i] = 0.0
    for i = 2:nx
        a[i] = 12.0/(dx*dx)-2.0/(α*dt)
        b[i] = -24.0/(dx*dx) - 20.0/(α*dt)
        c[i] = 12.0/(dx*dx)-2.0/(α*dt)
        r[i] = -2.0/(α*dt)*(u_n[k-1,i+1] + 10.0*u_n[k-1,i] + u_n[k-1,i-1]) -
                12.0/(dx*dx)*(u_n[k-1,i+1] - 2.0*u_n[k-1,i] + u_n[k-1,i-1])
    end
    i = nx+1
    a[i] = 0.0
    b[i] = 1.0
    c[i] = 0.0
    r[i] = 0.0
    u_n[k,:] = tdma(a,b,c,r,s,e)
end

u_error = zeros(nx+1)
u_error = u_n[nt+1,:] - u_e
rms_error = compute_l2norm(nx, u_error)

p1 = plot(x,u_e,lw = 4,xlabel="X", color = :red, ylabel = "U", xlims=(minimum(x),maximum(x)),
     grid=(:none), label = "Exact")

p2 = plot(x,u_n[nt+1,:],lw = 4,xlabel="X", color = :blue, ylabel = "U", xlims=(minimum(x),maximum(x)),
     grid=(:none), label = "t=1")

plot(p1, p2, size = (1000, 400))
