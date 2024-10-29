using Symbolics
using Roots
using DifferentialEquations
using Interpolations
using DataFrames
using CSV

# Parameters
mass = 1
spin = 3/10 * mass
Ei = 0.953
Lzi = 3.422
ri = 11.295
θi = π/2

# Variables for calculations
@variables M a
@variables τ
@variables t(τ) r(τ) θ(τ) ϕ(τ)
@variables ut ur uθ uϕ
ut = Symbolics.derivative(t, τ)
ur = Symbolics.derivative(r, τ)
uθ = Symbolics.derivative(θ, τ)
uϕ = Symbolics.derivative(ϕ, τ)

coords = [t, r, θ, ϕ] # Coordinates
fourvel = [ut, ur, uθ, uϕ] # 4-velocity
n = length(coords)

# Metric
Σ = r^2 + a^2 * cos(θ)^2
Δ = r^2 - 2 * M * r + a^2
tt = -1 + (2 * M * r) / Σ
rr = Σ / Δ
θθ = Σ
ϕϕ = sin(θ)^2 * (r^2 + a^2 + (2 * M * a^2 * r * sin(θ)^2) / Σ)
tϕ = -(2 * M * a * r * sin(θ)^2) / Σ
metric = [tt 0 0 tϕ; 0 rr 0 0; 0 0 θθ 0; tϕ 0 0 ϕϕ]
inversemetric = inv(metric)

# Christoffel symbols
Γ = [(1/2) * sum([inversemetric[i, s] * 
                (Symbolics.derivative(metric[s,j], coords[k]) + 
                    Symbolics.derivative(metric[s,k], coords[j]) - 
                    Symbolics.derivative(metric[j,k], coords[s])) 
                for s = 1:n]) 
        for i = 1:n, j = 1:n, k = 1:n]

# Geodesic equations
geodesic = [-sum([Γ[i,j,k] * Symbolics.derivative(coords[j], τ) * 
        Symbolics.derivative(coords[k], τ) for j = 1:n, k = 1:n]) 
    for i = 1:n]

# Returns the geodesic equations as numerical functions
function calculategeodesiceqns(mass, spin)
    geodesicfuncs = [build_function(substitute(geodesic[i], 
                                                Dict([M => mass, a => spin])), 
                                    ut, ur, uθ, uϕ, t, r, θ, ϕ) 
                    for i = 1:n]
    geodesicnums = [eval(geodesicfuncs[i]) for i = 1:n]

    return geodesicnums
end

function calculate4velocity(mass, spin, uti, uri, uϕi, ri, θi)
    # Solve for θ-component using normalization condition
    @variables uθi
    normalizationexpr = substitute(tt*uti^2 + rr*uri^2 + θθ*uθi^2 + ϕϕ*uϕi^2 + 2*tϕ*uti*uϕi + 1, 
                                    Dict([M => mass, a => spin, r => ri, θ => θi]))
    normalizationfunc(uθ) = substitute(normalizationexpr, Dict(uθi => uθ))
    uθi = find_zeros(normalizationfunc, -1e16, 1e16)[1]

    # Need to convert Num to Float64
    return [Float64(Symbolics.value(v)) for v in [uti, uri, uθi, uϕi]]
end

# Contravariant components of 4-velocity
uti = substitute((Ei*ϕϕ + Lzi*tϕ)/(tϕ^2 - tt*ϕϕ), 
                    Dict([M => mass, a => spin, r => ri, θ => θi]))
uri = 0
uϕi = substitute(-(Ei*tϕ + Lzi*tt)/(tϕ^2 - tt*ϕϕ), 
                    Dict([M => mass, a => spin, r => ri, θ => θi]))

# Initial conditions
initialpos = [0, ri, θi, 0]
initialvel = calculate4velocity(mass, spin, uti, uri, uϕi, ri, θi)                  

# Function to solve right-hand side of the geodesic equations
geodesics = calculategeodesiceqns(mass, spin)
function geodesiceqns!(ddu, du, u, p, τ)

    for i in 1:n
        ddu[i] = geodesics[i](du[1], du[2], du[3], du[4], u[1], u[2], u[3], u[4])
    end
    
end

geotspan = 10000
geoprob = SecondOrderODEProblem(geodesiceqns!, initialvel, initialpos, geotspan)
geodesicsoln = solve(geoprob, abstol=1e-10, reltol=1e-10)

# Solution
τsoln = geodesicsoln.t
tsoln = geodesicsoln[5, :]
rsoln = geodesicsoln[6, :]
θsoln = geodesicsoln[7, :]
ϕsoln = geodesicsoln[8, :]
utsoln = geodesicsoln[1, :]
ursoln = geodesicsoln[2, :]
uθsoln = geodesicsoln[3, :]
uϕsoln = geodesicsoln[4, :]
numpoints = length(τsoln)

# Interpolate and save solution
# followed example from https://github.com/JuliaMath/Interpolations.jl/

# Interpolation functions
interpolatet = linear_interpolation(τsoln, tsoln)
interpolater = linear_interpolation(τsoln, rsoln)
interpolateθ = linear_interpolation(τsoln, θsoln)
interpolateϕ = linear_interpolation(τsoln, ϕsoln)
interpolateut = linear_interpolation(τsoln, utsoln)
interpolateur = linear_interpolation(τsoln, ursoln)
interpolateuθ = linear_interpolation(τsoln, uθsoln)
interpolateuϕ = linear_interpolation(τsoln, uϕsoln)

τitp = range(τsoln[1], τsoln[numpoints], numpoints)
titp = interpolatet.(τitp) 
ritp = interpolater.(τitp) 
θitp = interpolateθ.(τitp) 
ϕitp = interpolateϕ.(τitp) 
utitp = interpolateut.(τitp) 
uritp = interpolateur.(τitp) 
uθitp = interpolateuθ.(τitp) 
uϕitp = interpolateuϕ.(τitp) 

df = DataFrame(τ = τitp, 
               t = titp, r = ritp, θ = θitp, ϕ = ϕitp,
               ut = utitp, ur = uritp, uθ = uθitp, uϕ = uϕitp)
CSV.write("geodesic-orbit.csv", df)