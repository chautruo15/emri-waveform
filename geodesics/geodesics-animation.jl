using DelimitedFiles
using Plots

geodesicorbit = DelimitedFiles.readdlm("geodesic-orbit.csv", ',', header=true)
orbit = geodesicorbit[1]

τ = orbit[:,1]
t = orbit[:,2]
r = orbit[:,3]
θ = orbit[:,4]
ϕ = orbit[:,5]
ut = orbit[:,6]
ur = orbit[:,7]
uθ = orbit[:,8]
uϕ = orbit[:,9]

n = length(τ)

# Convert from spherical to Cartesian coordinates
x = []
y = []
z = []
for i in eachindex(τ)
    push!(x, r[i] * sin(θ[i]) * cos(ϕ[i]))
    push!(y, r[i] * sin(θ[i]) * sin(ϕ[i]))
    push!(z, r[i] * cos(θ[i]))
end

# # Animation
# anim = @animate for i in eachindex(x)
#     plot3d(x[1:i], y[1:i], z[1:i], 
#         linewidth = 2,
#         xlabel = "x",
#         ylabel = "y",
#         zlabel = "z",
#         xlabelfontsize = 15,
#         ylabelfontsize = 15,
#         zlabelfontsize = 15,
#         xtickfontsize = 13,
#         ytickfontsize = 13,
#         ztickfontsize = 13,
#         xlim = (-13, 13),
#         ylim = (-13, 13),
#         zlim = (-3.5, 3.5),
#         size = (1000, 1000),
#         legend = false)
#     scatter!([x[i]], [y[i]], [z[i]], 
#         color=:red, 
#         markerstrokewidth=0, 
#         markersize=8)
# end

# gif(anim, "geodesics.gif", fps=50)
# gif(anim, "geodesics.mp4", fps=50)

# Animation
anim = @animate for i in eachindex(x)
    plot(x[1:i], y[1:i], 
        linewidth = 2,
        xlabel = "x",
        ylabel = "y",
        xlabelfontsize = 15,
        ylabelfontsize = 15,
        xtickfontsize = 13,
        ytickfontsize = 13,
        xlim = (-13, 13),
        ylim = (-13, 13),
        size = (1000, 1000),
        legend = false)
    scatter!([x[i]], [y[i]], 
        color=:red, 
        markerstrokewidth=0, 
        markersize=8)
end

gif(anim, "geodesics-2D.gif", fps=50)
gif(anim, "geodesics-2D.mp4", fps=50)