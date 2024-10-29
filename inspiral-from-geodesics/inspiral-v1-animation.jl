using DelimitedFiles
using Plots

# Read in orbit data
inspiralorbit = DelimitedFiles.readdlm("inspiral-orbit-v1.csv", ',', header=true)
orbit = inspiralorbit[1]

tins = orbit[:,1]
rins = orbit[:,2]
θins = orbit[:,3]
ϕins = orbit[:,4]
utins = orbit[:,5]
urins = orbit[:,6]
uθins = orbit[:,7]
uϕins = orbit[:,8]
utdotins = orbit[:,9]
urdotins = orbit[:,10]
uθdotins = orbit[:,11]
uϕdotins = orbit[:,12]

n = length(tins)

# Interpolate
interpolater = linear_interpolation(tins, rins)
interpolateθ = linear_interpolation(tins, θins)
interpolateϕ = linear_interpolation(tins, ϕins)
interpolateut = linear_interpolation(tins, utins)
interpolateur = linear_interpolation(tins, urins)
interpolateuθ = linear_interpolation(tins, uθins)
interpolateuϕ = linear_interpolation(tins, uϕins)
interpolateutdot = linear_interpolation(tins, utdotins)
interpolateurdot = linear_interpolation(tins, urdotins)
interpolateuθdot = linear_interpolation(tins, uθdotins)
interpolateuϕdot = linear_interpolation(tins, uϕdotins)

titp = range(tins[1], tins[n], 3000)
ritp = interpolater.(titp) 
θitp = interpolateθ.(titp) 
ϕitp = interpolateϕ.(titp) 
utitp = interpolateut.(titp) 
uritp = interpolateur.(titp) 
uθitp = interpolateuθ.(titp) 
uϕitp = interpolateuϕ.(titp) 
utdotitp = interpolateutdot.(titp) 
urdotitp = interpolateurdot.(titp) 
uθdotitp = interpolateuθdot.(titp) 
uϕdotitp = interpolateuϕdot.(titp) 

# Convert from spherical to Cartesian coordinates
x = [ritp[i] * sin(θitp[i]) * cos(ϕitp[i]) for i in eachindex(ritp)]
y = [ritp[i] * sin(θitp[i]) * sin(ϕitp[i]) for i in eachindex(ritp)]
z = [ritp[i] * cos(θitp[i]) for i in eachindex(ritp)]

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

# gif(anim, "inspiral-v1.gif", fps=50)
# gif(anim, "inspiral-v1.mp4", fps=50)

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

gif(anim, "inspiral-v1-2D.gif", fps=50)
gif(anim, "inspiral-v1-2D.mp4", fps=50)