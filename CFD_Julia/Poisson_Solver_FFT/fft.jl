clearconsole()

using CPUTime
using Printf
using Plots
using FFTW


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

nx = 64
ny = 64

dx = 2.0*pi/(nx)
dy = 2.0*pi/(ny)

x = 0:dx:2.0*pi
y = 0:dy:2.0*pi

ue = Array{Float64}(undef,nx+1,ny+1)
f = Array{Float64}(undef,nx+1,ny+1)
un = Array{Float64}(undef,nx+1,ny+1)

# test case Ue = sin(3x) + cos(2y); f = -9sin(3x) - 4cos(2y)
for j = 1:ny+1
    for i = 1:nx+1
        ue[i,j] = sin(3.0*x[i]) + cos(2.0*y[j])
        f[i,j] = -9.0*sin(3.0*x[i]) - 4.0*cos(2.0*y[j])
    end
end

un[1:nx,1:ny] = fps(nx,ny,dx,dy,f)

# extend the boundary of the domain
un[:,ny+1] = un[:,1]
un[nx+1,:] = un[1,:]

p1 = contour(x, y, transpose(ue), fill=true)
p2 = contour(x, y, transpose(un), fill=true)
p3 = plot(p1,p2, size = (1000, 400))
savefig(p3,"contour.pdf")
