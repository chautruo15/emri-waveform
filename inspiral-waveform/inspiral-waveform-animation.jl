using DelimitedFiles
using Interpolations
using Plots
using LaTeXStrings

# Read in data
inspiralorbit = DelimitedFiles.readdlm("inspiral-e0.0-p20.0-iota1.57.csv", ',', header=true)
orbit = inspiralorbit[1]
t = orbit[:,1]
r = orbit[:,2]
θ = orbit[:,3]
ϕ = orbit[:,4]

waveform = DelimitedFiles.readdlm("waveform-e0.0-p20.0-iota1.57.csv", ',', header=true)
wave = waveform[1]
hplusdata = wave[:,1]
hcrossdata = wave[:,2]

n = length(t)

# Interpolate

m = 200000 # New number of points

# Inspiral orbit
itpr = linear_interpolation(t, r)
itpθ = linear_interpolation(t, θ)
itpϕ = linear_interpolation(t, ϕ)

titp = range(t[1], t[n], m)
ritp = itpr.(titp) 
θitp = itpθ.(titp) 
ϕitp = itpϕ.(titp) 

# Waveform
itphplus = linear_interpolation(t, hplusdata)
itphcross = linear_interpolation(t, hcrossdata)

hplusitp = itphplus.(titp)
hcrossitp = itphcross.(titp)

# Convert to Cartesian
xitp = [ritp[i] * sin(θitp[i]) * cos(ϕitp[i]) for i in 1:m]
yitp = [ritp[i] * sin(θitp[i]) * sin(ϕitp[i]) for i in 1:m]
zitp = [ritp[i] * cos(θitp[i]) for i in 1:m]

start = 2500
numpoints = 12000
stop = start + numpoints - 1

# Two-panel animation
anim = @animate for i in (start+1):stop
    if (i <= 250)
        p1 = plot3d(xitp[start:i], yitp[start:i], zitp[start:i], 
            linewidth = 1.75,
            xlabel = "x",
            ylabel = "y",
            zlabel = "z",
            xlabelfontsize = 15,
            ylabelfontsize = 15,
            zlabelfontsize = 15,
            xtickfontsize = 13,
            ytickfontsize = 13,
            ztickfontsize = 13,
            xlim = (-25, 25),
            ylim = (-25, 25),
            zlim = (-25, 25),
            seriesalpha = range(0, 1, length = length(xitp[start:i])),
            size = (800, 800),
            legend = false)

        scatter!(p1, [xitp[i]], [yitp[i]], [zitp[i]], 
            color=:red, 
            markerstrokewidth=0, 
            markersize=8)
    else
        p1 = plot3d(xitp[i - 249:i], yitp[i - 249:i], zitp[i - 249:i], 
            linewidth = 1.75,
            xlabel = "x",
            ylabel = "y",
            zlabel = "z",
            xlabelfontsize = 15,
            ylabelfontsize = 15,
            zlabelfontsize = 15,
            xtickfontsize = 13,
            ytickfontsize = 13,
            ztickfontsize = 13,
            xlim = (-25, 25),
            ylim = (-25, 25),
            zlim = (-25, 25),
            seriesalpha = range(0, 1, length = 250),
            size = (800, 800),
            legend = false)

        scatter!(p1, [xitp[i]], [yitp[i]], [zitp[i]], 
            color=:red, 
            markerstrokewidth=0, 
            markersize=8)
    end
    
    p2 = plot(titp[start:i], hplusitp[start:i], 
        linewidth = 1.75,
        xlabel = "t",
        ylabel = L"h_+",
        xaxis = nothing,
        yaxis = nothing,
        size = (1000, 200),
        legend = false)
    
    scatter!(p2, [titp[i]], [hplusitp[i]],
        color=:red, 
        markerstrokewidth=0, 
        markersize=4)
    
    plot(p1, p2, layout = grid(2, 1, heights=[5/6, 1/6]), size = (1000, 1000))
end

gif(anim, "inspiral-waveform-animation-e0.0-p20.0-iota1.57-v2.gif", fps=200)
gif(anim, "inspiral-waveform-animation-e0.0-p20.0-iota1.57-v2.mp4", fps=200)