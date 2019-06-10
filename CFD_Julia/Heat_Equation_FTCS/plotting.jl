clearconsole()

using CSV
using Plots
font = Plots.font("Times New Roman", 20)
pyplot(guidefont=font, xtickfont=font, ytickfont=font, legendfont=font)

final_field = CSV.read("field_final.csv")#, datarow = 2, type=Float64)

x = convert(Array,final_field[:,1])

u_e = convert(Array,final_field[:,2])
u_n = convert(Array,final_field[:,3])
u_error = convert(Array,final_field[:,4])

u = Array{Float64}(undef, length(u_e), 2)
u[:,1] = u_e
u[:,2] = u_n

# plotting the exact and numerical solution
styles = [:solid, :dot]
styles = reshape(styles, 1, length(styles))
color=[:red :blue]

p1 = plot(x,u, line = (6,styles), xlabel="\$X\$", ylabel = "\$U\$",  xlims=(minimum(x),maximum(x)),
     grid=(:none), linestyle=:dot, label=["Exact" "Numerical"], color=color)
     #legendfontsize=15, xtickfontsize=15, ytickfontsize=15, guidefontsize=15)

plot(p1, size = (600, 400))
savefig(p1,"ftcs.pdf")
