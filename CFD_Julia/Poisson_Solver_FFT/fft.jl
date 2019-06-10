clearconsole()

using CPUTime
using Printf
using Plots
using FFTW

function compute_l2norm(nx, ny, r)
    rms = 0.0
    # println(residual)
    for j = 2:ny for i = 2:nx
        rms = rms + r[i,j]^2
    end end
    # println(rms)
    rms = sqrt(rms/((nx-1)*(ny-1)))
    return rms
end

function fps(nx,ny,dx,dy,f)
    eps = 1.0e-6

    kx = Array{Float64}(undef,nx)
    ky = Array{Float64}(undef,ny)

    data = Array{Complex{Float64}}(undef,nx,ny)
    data1 = Array{Complex{Float64}}(undef,nx,ny)
    e = Array{Complex{Float64}}(undef,nx,ny)

    u = Array{Complex{Float64}}(undef,nx,ny)

    aa = -2.0/(dx*dx) - 2.0/(dy*dy)
    bb = 2.0/(dx*dx)
    cc = 2.0/(dy*dy)

    hx = 2.0*pi/nx
    hy = 2.0*pi/ny

    for i = 1:nx
        kx[i] = hx*(i-1.0)
    end

    for i = 1:ny
        ky[i] = hy*(i-1.0)
    end

    kx[1] = eps
    ky[1] = eps

    for i = 1:nx
        for j = 1:ny
            data[i,j] = complex(f[i,j],0.0)
        end
    end

    e = fft(data)
    e[1,1] = 0.0
    for i = 1:nx
        for j = 1:ny
            data1[i,j] = e[i,j]/(aa + bb*cos(kx[i]) + cc*cos(ky[j]))
        end
    end

    u = real(ifft(data1))

    return u
end

nx = 128
ny = 128

# dx = 2.0*pi/(nx)
# dy = 2.0*pi/(ny)
# x = 0:dx:2.0*pi
# y = 0:dy:2.0*pi

x_l = 0.0
x_r = 1.0
y_b = 0.0
y_t = 1.0

x = Array{Float64}(undef, nx+1)
y = Array{Float64}(undef, ny+1)
ue = Array{Float64}(undef,nx+1,ny+1)
f = Array{Float64}(undef,nx+1,ny+1)
un = Array{Float64}(undef,nx+1,ny+1)

for i = 1:nx+1
    x[i] = x_l + dx*(i-1)
end
for i = 1:ny+1
    y[i] = y_b + dy*(i-1)
end

c1 = (1.0/16.0)^2
c2 = -2.0*pi*pi

# test case Ue = sin(3x) + cos(2y); f = -9sin(3x) - 4cos(2y)
for j = 1:ny+1
    for i = 1:nx+1
        # ue[i,j] = sin(3.0*x[i]) + cos(2.0*y[j])
        # f[i,j] = -9.0*sin(3.0*x[i]) -4.0*cos(2.0*y[j])
        # un[i,j] = 0.0
        ue[i,j] = sin(pi*x[i])*sin(pi*y[j]) +
                  c1*sin(16.0*pi*x[i])*sin(16.0*pi*y[j])

        f[i,j] = c2*sin(pi*x[i])*sin(pi*y[j]) +
                 c2*sin(16.0*pi*x[i])*sin(16.0*pi*y[j])

        un[i,j] = 0.0
    end
end

# Dirichlet boundary condition
un[:,1] = ue[:,1]
un[:, ny+1] = ue[:, ny+1]

un[1,:] = ue[1,:]
un[nx+1,:] = ue[nx+1,:]

val, t, bytes, gctime, memallocs = @timed begin

un[1:nx,1:ny] = fps(nx,ny,dx,dy,f)

end

# Dirichlet boundary condition
un[:,1] = ue[:,1]
un[:, ny+1] = ue[:, ny+1]

un[1,:] = ue[1,:]
un[nx+1,:] = ue[nx+1,:]

uerror = zeros(nx+1, ny+1)
rms_error = 0.0

uerror = un - ue

rms_error = compute_l2norm(nx, ny, uerror)
max_error = maximum(abs.(uerror))

println("Error details:");
println("L-2 Norm = ", rms_error);
println("Maximum Norm = ", max_error);
print("CPU Time = ", t);

p1 = contour(x, y, transpose(ue), fill=true)
p2 = contour(x, y, transpose(un), fill=true)
p3 = plot(p1,p2, size = (1000, 400))
savefig(p3,"contour.pdf")
