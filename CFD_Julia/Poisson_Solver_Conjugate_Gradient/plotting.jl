clearconsole()

using CSV
using Plots
font = Plots.font("Times New Roman", 18)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

nx = 128
ny = 128
residual_hist = CSV.read("residual.csv")#, datarow = 3, type=Float64)
iter_hist =  residual_hist[:,1]
res_hist = convert(Matrix, residual_hist[:,2:3])

p = plot(iter_hist,res_hist,lw = 2,yscale = :log10, xlabel="Iteration count",
               ylabel = "L2 Norm", xlims=(0,maximum(iter_hist)),grid=(:none),
               label=["rms" "rms/rms\$_0\$"])

savefig(p,"myplot.pdf")

init_field = CSV.read("field_initial.csv")#, type=Float64)
final_field = CSV.read("field_final.csv")#, datarow = 3, type=Float64)

x = convert(Array,init_field[:,1])
y = convert(Array,init_field[:,2])

u_e = convert(Array,init_field[:,5])
u_e = reshape(u_e, (nx+1, ny+1))

u_n = convert(Array,final_field[:,4])
u_n = reshape(u_n, (nx+1, ny+1))

xx = x[1:nx+1]
yy = reshape(y, (nx+1, ny+1))[1,:]

XX = repeat(xx, 1, length(yy))
XX = convert(Matrix,transpose(XX))
YY = repeat(yy, 1, length(xx))

p1 = contour(xx, yy, u_e, fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Exact")
p2 = contour(xx, yy, u_n, fill=true,xlabel="\$X\$", ylabel="\$Y\$", title="Numerical")
p3 = plot(p1,p2, size = (1300, 600))
savefig(p3,"contour.pdf")
