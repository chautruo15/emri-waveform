using DelimitedFiles
using DataFrames
using CSV

# Read in orbit data
inspiralorbit = DelimitedFiles.readdlm("inspiral.csv", ',', header=true)
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

# Convert r, θ, ϕ to x, y, z
function calculatexyz(rsoln, θsoln, ϕsoln)
    m = length(rsoln)

    x = [rsoln[i] * sin(θsoln[i]) * cos(ϕsoln[i]) for i in 1:m]
    y = [rsoln[i] * sin(θsoln[i]) * sin(ϕsoln[i]) for i in 1:m]
    z = [rsoln[i] * cos(θsoln[i]) for i in 1:m]

    return x, y, z
end

# First derivative of x, y, z wrt τ
function calculatexyzτdot(rsoln, θsoln, ϕsoln, 
                        ursoln, uθsoln, uϕsoln)
    m = length(rsoln)

    xτdot = [cos(ϕsoln[i]) * ursoln[i] * sin(θsoln[i]) + 
            uθsoln[i] * cos(ϕsoln[i]) * cos(θsoln[i]) * rsoln[i] -
            sin(θsoln[i]) * uϕsoln[i] * sin(ϕsoln[i]) * rsoln[i]
        for i in 1:m]

    yτdot = [ursoln[i] * sin(θsoln[i]) * sin(ϕsoln[i]) + 
            uθsoln[i] * cos(θsoln[i]) * sin(ϕsoln[i]) * rsoln[i] +
            cos(ϕsoln[i]) * sin(θsoln[i]) * uϕsoln[i] * rsoln[i]
        for i in 1:m]

    zτdot = [ursoln[i] * cos(θsoln[i]) - uθsoln[i] * sin(θsoln[i]) * rsoln[i]
        for i in 1:m]

    return xτdot, yτdot, zτdot
end

# second derivative of x, y, z wrt τ
function calculatexyzτdotdot(rsoln, θsoln, ϕsoln, 
                            ursoln, uθsoln, uϕsoln,
                            urdotsoln, uθdotsoln, uϕdotsoln)
    m = length(rsoln)

    xτdotdot = [cos(ϕsoln[i]) * urdotsoln[i] * sin(θsoln[i]) + 
                2 * uθsoln[i] * cos(ϕsoln[i]) * ursoln[i] * cos(θsoln[i]) +
                cos(ϕsoln[i]) * uθdotsoln[i] * cos(θsoln[i]) * rsoln[i] -
                ursoln[i] * sin(θsoln[i]) * uϕsoln[i] * sin(ϕsoln[i]) -
                sin(θsoln[i]) * uϕdotsoln[i] * sin(ϕsoln[i]) * rsoln[i] -
                uθsoln[i]^2 * cos(ϕsoln[i]) * sin(θsoln[i]) * rsoln[i] -
                2 * uθsoln[i] * uϕsoln[i] * cos(θsoln[i]) * sin(ϕsoln[i]) * rsoln[i] -
                uϕsoln[i]^2 * cos(ϕsoln[i]) * sin(θsoln[i]) * rsoln[i]
        for i in 1:m]

    yτdotdot = [urdotsoln[i] * sin(θsoln[i]) * sin(ϕsoln[i]) +
                2 * uθsoln[i] * ursoln[i] * cos(θsoln[i]) * sin(ϕsoln[i]) +
                2 * cos(ϕsoln[i]) * ursoln[i] * sin(θsoln[i]) * uϕsoln[i] +
                cos(ϕsoln[i]) * sin(θsoln[i]) * uϕdotsoln[i] * rsoln[i] +
                uθdotsoln[i] * cos(θsoln[i]) * sin(ϕsoln[i]) * rsoln[i] -
                uθsoln[i]^2 * sin(θsoln[i]) * sin(ϕsoln[i]) * rsoln[i] +
                2 * uθsoln[i] * cos(ϕsoln[i]) * uϕsoln[i] * cos(θsoln[i]) * rsoln[i] -
                uϕsoln[i]^2 * sin(θsoln[i]) * sin(ϕsoln[i]) * rsoln[i]
        for i in 1:m]

    zτdotdot = [urdotsoln[i] * cos(θsoln[i]) -
                2 * uθsoln[i] * ursoln[i] * sin(θsoln[i]) -
                sin(θsoln[i]) * uθdotsoln[i] * rsoln[i] -
                uθsoln[i]^2 * cos(θsoln[i]) * rsoln[i]
        for i in 1:m]

    return xτdotdot, yτdotdot, zτdotdot
end

function calculatexyztdot(xτdot, yτdot, zτdot, utsoln)
    m = length(xτdot)

    xtdot = [xτdot[i] / utsoln[i] for i in 1:m]
    ytdot = [yτdot[i] / utsoln[i] for i in 1:m]
    ztdot = [zτdot[i] / utsoln[i] for i in 1:m]

    return xtdot, ytdot, ztdot
end

function calculatexyztdotdot(xτdotdot, yτdotdot, zτdotdot, utsoln)
    m = length(xτdotdot)

    xtdotdot = [xτdotdot[i] / utsoln[i] for i in 1:m]
    ytdotdot = [yτdotdot[i] / utsoln[i] for i in 1:m]
    ztdotdot = [zτdotdot[i] / utsoln[i] for i in 1:m]

    return xtdotdot, ytdotdot, ztdotdot
end

function convertpctonaturalunits(mass, distance)
    conversionfactor = 3.0857e16 # m/pc
    c = 2.998e8 # m/s 
    G = 6.674e-11 # m^3/kg/s^2
    scalemass = 1e6
    Msolar = 1.989e30 # kg/solar mass

    return (distance * conversionfactor * c^2) / (G * mass * scalemass * Msolar)
end

function calculatehbar(massratio, distance, xyz, xyztdot, xyztdotdot)
    Idotdot = [massratio * xyztdotdot[j] * xyz[k] + 
            massratio * xyz[j] * xyztdotdot[k] + 
            2 * massratio * xyztdot[j] * xyztdot[k]
        for j = 1:3, k = 1:3]

    hbar = (2/distance) * Idotdot 

    return hbar
end

function calculatehcomponents(Θ, Φ, hbar)
    hΘΘ = cos(Θ)^2 * (hbar[1,1]*cos(Φ)^2 + hbar[1,2]*sin(2*Φ) + hbar[2,2]*sin(Φ)^2) +
    hbar[3,3]*sin(Θ)^2 - sin(2*Θ) * (hbar[1,3]*cos(Φ) + hbar[2,3]*sin(Φ))
    hΘΦ = cos(Θ) * (-(1/2)*hbar[1,1]*sin(2*Φ) + hbar[1,2]*cos(2*Φ) + (1/2)*hbar[2,2]*sin(2*Φ)) +
        sin(Θ) * (hbar[1,3]*sin(Φ) - hbar[2,3]*cos(Φ))
    hΦΦ = hbar[1,1]*sin(Φ)^2 - hbar[1,2]*sin(2*Φ) + hbar[2,2]*cos(Φ)^2

    hPlus = hΘΘ - hΦΦ
    hCross = 2 * hΘΦ

    return hPlus, hCross
end

# Parameters
mass = 1    # in units of 1e6 solar masses
massratio = 1e-2
distance = convertpctonaturalunits(mass, 5e9) # 5e9 pc
bigtheta = π/4
bigphi = 0

# To save waveform components
hPluses = Vector{Float64}()
hCrosses = Vector{Float64}()

x, y, z = calculatexyz(rins, θins, ϕins)
xτdot, yτdot, zτdot = calculatexyzτdot(rins, θins, ϕins, urins, uθins, uϕins)
xτdotdot, yτdotdot, zτdotdot = calculatexyzτdotdot(rins, θins, ϕins, 
                                                urins, uθins, uϕins,
                                                urdotins, uθdotins, uϕdotins)

xtdot, ytdot, ztdot = calculatexyztdot(xτdot, yτdot, zτdot, utins)
xtdotdot, ytdotdot, ztdotdot = calculatexyztdotdot(xτdotdot, yτdotdot, zτdotdot, utins)

for i = 1:length(tins)
    xyzarr = [x[i], y[i], z[i]]
    xyzdotarr = [xtdot[i], ytdot[i], ztdot[i]]
    xyzdotdotarr = [xtdotdot[i], ytdotdot[i], ztdotdot[i]]

    hbar = calculatehbar(massratio, distance, xyzarr, xyzdotarr, xyzdotdotarr)
    hPlus, hCross = calculatehcomponents(bigtheta, bigphi, hbar)

    push!(hPluses, hPlus)
    push!(hCrosses, hCross)
end

# Save waveform
df = DataFrame(plus = hPluses, cross = hCrosses)
CSV.write("waveform.csv", df)