using Symbolics
using Roots
using DifferentialEquations
using DataFrames
using CSV
using Interpolations

# Parameters
mass = 1    # in units of 1e6 solar masses
spin = 3/10 * mass
littlemass = 1
massratio = 1e-2
e = 0
p = 20
ι = π/2 - 0.001
saveepιfluxes = false   # Specify whether to save orbital parameters

# To save inspiral data
tinspiral = Vector{Float64}()
rinspiral = Vector{Float64}()
θinspiral = Vector{Float64}()
ϕinspiral = Vector{Float64}()
utinspiral = Vector{Float64}()
urinspiral = Vector{Float64}()
uθinspiral = Vector{Float64}()
uϕinspiral = Vector{Float64}()
utdotinspiral = Vector{Float64}()
urdotinspiral = Vector{Float64}()
uθdotinspiral = Vector{Float64}()
uϕdotinspiral = Vector{Float64}()

if saveepιfluxes
    tfluxes = Vector{Float64}()
    einspiral = Vector{Float64}()
    pinspiral = Vector{Float64}()
    ιinspiral = Vector{Float64}()
    Efluxes = Vector{Float64}()
    Lfluxes = Vector{Float64}()
    Qfluxes = Vector{Float64}()
end

# Variables for calculations
@variables M a μ
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
Σ = r^2 + a^2*cos(θ)^2
Δ = r^2 - 2*M*r +a^2
tt = -1 + (2*M*r)/Σ
rr = Σ/Δ
θθ = Σ
ϕϕ = sin(θ)^2 * (r^2 + a^2 + (2*M*a^2*r*sin(θ)^2)/Σ)
tϕ = -(2*M*a*r*sin(θ)^2)/Σ
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

function calculatehorizon(mass, spin)
    return mass + sqrt(mass^2 - spin^2)
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

function calculatelastvelocityposition(soln)
    m = length(soln.t)

    utf = soln[1, :][m]
    urf = soln[2, :][m]
    uθf = soln[3, :][m]
    uϕf = soln[4, :][m]
    tf = soln[5, :][m]
    rf = soln[6, :][m]
    θf = soln[7, :][m]
    ϕf = soln[8, :][m]

    return [utf, urf, uθf, uϕf], [tf, rf, θf, ϕf]
end

function calculateepι(soln)
    rsoln = soln[6, :]
    θsoln = soln[7, :]

    ra = maximum(rsoln) # Apastron (max)
    rp = minimum(rsoln) # Periastron (min)
    
    e = (ra - rp) / (ra + rp) # Eccentricity
    p = (2*ra*rp) / (ra + rp) # Semilatus
    ι = π/2 - minimum(θsoln) # Inclination angle

    return e, p, ι
end

# Coefficients to solve for E, Lz, and Q from e, p, and ι
# For non-circular orbits
αI(mass, spin, rI, zm) = (rI^2 + spin^2) * (rI^2 + spin^2 * zm) + 2 * mass * rI * spin^2 * (1.0 - zm) # Eq. E4
βI(mass, spin, rI, zm) = -2.0 * mass * rI * spin    # Eq. E5
γI(mass, spin, rI, zm) = -(1.0 / (1 - zm)) * (rI^2 + spin^2 * zm - 2.0 * mass * rI) # Eq. E6
λI(mass, spin, rI, zm) = -(rI^2 + spin^2 * zm) * (rI^2 - 2.0 * mass * rI + spin^2) # Eq. E7

# For circular orbits
α_2(mass, spin, r0, zm) = 2.0 * r0 * (r0^2 + spin^2) - spin^2 * (r0 - mass) * (1.0 - zm)    # Eq. E8
β_2(mass, spin, r0, zm) = -spin * mass   # Eq. E9
γ_2(mass, spin, r0, zm) = -(r0 - mass) / (1.0 - zm)   # Eq. E10
λ_2(mass, spin, r0, zm) = -r0 * (r0^2 - 2.0 * mass * r0 + spin^2) - (r0 - mass) * (r0^2 + spin^2 * zm) # Eq. E11

commute(Πa, Πp, Ωa, Ωp) = Πa * Ωp - Πp * Ωa # Anti-symmetric product (Eq. E13)

# Returns conserved quantities
function calculateELQfromepι(mass, spin, e, p, ι)
    if (ι < π/2)
        signLz = 1
    else
        signLz = -1
    end

    θmin = π/2 - signLz * ι     # Eq. E1
    zm = cos(θmin)^2            # Eq. 92

    if (e == 0.0)
        r0 = p * mass
        α1 = αI(mass, spin, r0, zm)
        α2 = α_2(mass, spin, r0, zm)
        β1 = βI(mass, spin, r0, zm)
        β2 = β_2(mass, spin, r0, zm)
        γ1 = γI(mass, spin, r0, zm)
        γ2 = γ_2(mass, spin, r0, zm)
        λ1 = λI(mass, spin, r0, zm)
        λ2 = λ_2(mass, spin, r0, zm)
    else
        rp = (p * mass) / (1.0 + e)
        ra = (p * mass) / (1.0 - e)
        α1 = αI(mass, spin, ra, zm)
        α2 = αI(mass, spin, rp, zm)
        β1 = βI(mass, spin, ra, zm)
        β2 = βI(mass, spin, rp, zm)
        γ1 = γI(mass, spin, ra, zm)
        γ2 = γI(mass, spin, rp, zm)
        λ1 = λI(mass, spin, ra, zm)
        λ2 = λI(mass, spin, rp, zm)
    end

    # Coefficients of Eq. E12
    A = commute(α1, α2, γ1, γ2)^2 + 4 * commute(α1, α2, β1, β2) * commute(γ1, γ2, β1, β2)
    B = 2 * (commute(α1, α2, γ1, γ2) * commute(λ1, λ2, γ1, γ2) + 2.0 * commute(γ1, γ2, β1, β2) * commute(λ1, λ2, β1, β2))
    C = commute(λ1, λ2, γ1, γ2)^2

    if (signLz > 0)
        E = sqrt((-B - sqrt(B^2 - 4.0 * A * C)) / (2.0 * A))
        Lz = sqrt((commute(α1, α2, β1, β2) * E^2 + commute(λ1, λ2, β1, β2)) / commute(β1, β2, γ1, γ2))
    else
        E = sqrt((-B + sqrt(B^2 - 4.0 * A * C)) / (2.0 * A))
        Lz = -sqrt((commute(α1, α2, β1, β2) * E^2 + commute(λ1, λ2, β1, β2)) / commute(β1, β2, γ1, γ2))
    end

    if (θmin == 0.0)
        C = 0.0
    else    
        C = zm * (Lz^2 / (1.0 - zm) + spin^2 * (1.0 - E^2))
    end

    Q = C + (Lz - spin * E)^2  # Eq. 17

    return E, Lz, Q
end

function calculateELQfromvelpos(mass, spin, littlemass, finalvel, finalpos)
    utf = finalvel[1]
    uθf = finalvel[3]
    uϕf = finalvel[4]
    
    rf = finalpos[2]
    θf = finalpos[3]

    Ef = substitute(-tt*ut - tϕ*uϕ, Dict([M => mass, a => spin, 
                                            r => rf, θ => θf, 
                                            ut => utf, uϕ => uϕf]))
    Lzf = substitute(tϕ*ut + ϕϕ*uϕ, Dict([M => mass, a => spin, 
                                            r => rf, θ => θf, 
                                            ut => utf, uϕ => uϕf]))
    Qf = substitute((θθ*uθ)^2 + cos(θ)^2 * (a^2 * ((μ/M)^2 - Ef^2) + (Lzf/sin(θ))^2), 
                    Dict([M => mass, μ => littlemass, a => spin, 
                            r => rf, θ => θf, 
                            uθ => uθf]))

    # Need to convert Num to Float64                        
    return Float64(Symbolics.value(Ef)), 
        Float64(Symbolics.value(Lzf)), 
        Float64(Symbolics.value(Qf))
end

function calculateELQfluxes(mass, spin, massratio, e, p, ι, Q)
    g1 = 1 + (73/24)e^2 + (37/96)e^4
    g2 = (73/12) + (823/24)e^2 + (949/32)e^4 + (491/192)e^6
    g3 = (1247/336) + (9181/672)e^2
    g4 = 4 + (1375/48)e^2
    g5 = (44711/9072) + (172157/2592)e^2
    g6 = (33/16) + (359/32)e^2
    g7 = (8191/672) + (44531/336)e^2
    g8 = (3749/336) - (5143/168)e^2
    g9 = 1 + (7/8)e^2
    g10 = (61/12) + (119/8)e^2 + (183/32)e^4
    g11 = (1247/336) + (425/336)e^2
    g12 = 4 + (97/8)e^2
    g13 = (44711/9072) + (302893/6048)e^2
    g14 = (33/16) + (95/16)e^2
    g15 = (8191/672) + (48361/1344)e^2
    g16 = (417/56) - (37241/672)e^2

    g10a = (61/24) + (63/8)e^2 + (95/64)e^4
    g10b = (61/8) + (91/4)e^2 + (461/64)e^4

    # Eq. 44
    Eflux = -(32/5) * massratio/(p^5) * (1-e^2)^(3/2) * 
        (g1 - (a/M) * (M/p)^(3/2) * g2 * cos(ι) - (M/p) * g3 +
        π * (M/p)^(3/2) * g4 - (M/p)^2 * g5 + (a/M)^2 * (M/p)^2 * g6 -
        (527/96) * (a/M)^2 * (M/p)^2 * sin(ι)^2)
    # Eflux = -(32/5) * (μ/M)^2 * (M/p)^5 * (1-e^2)^(3/2) * 
    #     (g1 - (a/M) * (M/p)^(3/2) * g2 * cos(ι) - (M/p) * g3 +
    #     π * (M/p)^(3/2) * g4 - (M/p)^2 * g5 + (a/M)^2 * (M/p)^2 * g6 -
    #     (527/96) * (a/M)^2 * (M/p)^2 * sin(ι)^2)
    Eflux = substitute(Eflux, Dict([M => mass, a => spin]))

    # Eq. 45
    Lflux = -(32/5) * massratio/(M^(1/2) * p^(7/2)) * (1-e^2)^(3/2) *
        (g9 * cos(ι) + (a/M) * (M/p)^(3/2) * (g10a - cos(ι)^2 * g10b) -
        (M/p) * g11 * cos(ι) + π * (M/p)^(3/2) * g12 * cos(ι) -
        (M/p)^2 * g13 * cos(ι) + (a/M)^2 * (M/p)^2 * cos(ι) * 
        (g14 - (45/8) * sin(ι)^2))
    # Lflux = -(32/5) * (μ^2/M) * (M/p)^(7/2) * (1-e^2)^(3/2) *
    #     (g9 * cos(ι) + (a/M) * (M/p)^(3/2) * (g10a - cos(ι)^2 * g10b) -
    #     (M/p) * g11 * cos(ι) + π * (M/p)^(3/2) * g12 * cos(ι) -
    #     (M/p)^2 * g13 * cos(ι) + (a/M)^2 * (M/p)^2 * cos(ι) * 
    #     (g14 - (45/8) * sin(ι)^2))
    Lflux = substitute(Lflux, Dict([M => mass, a => spin]))

    # Eq. 56
    Qflux = -(64/5) * massratio/(M^(1/2) * p^(7/2)) * sqrt(Q) * sin(ι) * (1-e^2)^(3/2) * 
        (g9 - (a/M) * (M/p)^(3/2) * cos(ι) * g10b - (M/p) * g11 +
        π * (M/p)^(3/2) * g12 - (M/p)^2 * g13 + (a/M)^2 * (M/p)^2 *
        (g14 - (45/8) * sin(ι)^2))
    # Qflux = -(64/5) * (μ^2/M) * (M/p)^(7/2) * sqrt(Q) * sin(ι) * (1-e^2)^(3/2) * 
    #     (g9 - (a/M) * (M/p)^(3/2) * cos(ι) * g10b - (M/p) * g11 +
    #     π * (M/p)^(3/2) * g12 - (M/p)^2 * g13 + (a/M)^2 * (M/p)^2 *
    #     (g14 - (45/8) * sin(ι)^2))
    Qflux = substitute(Qflux, Dict([M => mass, a => spin]))

    # Need to convert Num to Float64
    return Float64(Symbolics.value(Eflux)), 
        Float64(Symbolics.value(Lflux)), 
        Float64(Symbolics.value(Qflux))
end

function calculateradiationforce(mass, spin, massratio, finalvel, finalpos, E, Lz, Eflux, Lflux, Qflux)
    utf = finalvel[1]
    urf = finalvel[2]
    uθf = finalvel[3]
    uϕf = finalvel[4]
    
    tf = finalpos[1]
    rf = finalpos[2]
    θf = finalpos[3]
    ϕf = finalpos[4]
    
    # Ft and Fϕ
    if abs(substitute(tϕ*tϕ - tt*ϕϕ, Dict([M => mass, a => spin,
        r => rf, θ => θf]))) < 1e-8 
        Ft = 0.0
        Fϕ = 0.0
    else
        A = [-tt -tϕ; tϕ ϕϕ]
        b = [Eflux*ut, Lflux*ut]
        Ft = (A\b)[1]
        Fϕ = (A\b)[2]
    end

    if abs(utf) < 1e-3
        Fθ = 0.0
    else
        @variables Fθ # solve for this
        Q̇eqn = 2*θθ^2*uθ*Fθ + 2*cos(θ)^2*a^2*E*Eflux + (2*cos(θ)^2*Lz*Lflux)/sin(θ)^2 ~ Qflux*ut
        Fθ = symbolic_linear_solve(Q̇eqn, Fθ)
    end

    if abs(urf) < (1e-1)*massratio || abs(urf) < 1e-16
        Fr = 0.0
    else
        @variables Fr # Solve for this
        eqn = tt*ut*Ft + tϕ*uϕ*Ft + rr*ur*Fr + θθ*uθ*Fθ + tϕ*ut*Fϕ + ϕϕ*uϕ*Fϕ ~ 0
        Fr = symbolic_linear_solve(eqn, Fr)
    end

    # Radiation force
    F = substitute([Ft, Fr, Fθ, Fϕ], Dict([M => mass, a => spin,
                                            t => tf, r => rf, θ => θf, ϕ => ϕf,
                                            ut => utf, ur => urf, uθ => uθf, uϕ => uϕf]))
    
    return F
end

# Returns the inspiral equations as numerical functions
function calculateinspiraleqns(mass, spin, radiationforce)
    inspiral = [substitute(geodesic[i], Dict([M => mass, a => spin])) + 
                radiationforce[i] for i = 1:n]
    inspiralfuncs = [build_function(inspiral[i], ut, ur, uθ, uϕ, t, r, θ, ϕ) for i = 1:n]
    inspiralnums = [eval(inspiralfuncs[i]) for i = 1:n]
    
    return inspiralnums
end

# Function to solve right-hand side of inspiral equations
function inspiraleqns!(ddu, du, u, p, τ)
    inspirals = p[1]
    n = length(inspirals)
    
    for i in 1:n
        ddu[i] = inspirals[i](du[1], du[2], du[3], du[4], u[1], u[2], u[3], u[4])
    end
end

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

ri = p * mass
E, Lz, Q = calculateELQfromepι(mass, spin, e, p, ι)
horizon = calculatehorizon(mass, spin)

# Contravariant components of 4-velocity
uti = substitute((E*ϕϕ + Lz*tϕ)/(tϕ^2 - tt*ϕϕ), 
                    Dict([M => mass, a => spin, r => ri, θ => ι]))
uri = 0
uϕi = substitute(-(E*tϕ + Lz*tt)/(tϕ^2 - tt*ϕϕ), 
                    Dict([M => mass, a => spin, r => ri, θ => ι]))

# Initial conditions
initialpos = [0, ri, ι, 0]
initialvel = calculate4velocity(mass, spin, uti, uri, uϕi, ri, ι)  

initialvelnew = initialvel 
initialposnew = initialpos
while (round(initialposnew[2], digits=5) > round(1.05*horizon, digits=5)) 
    Eflux, Lflux, Qflux = calculateELQfluxes(mass, spin, massratio, e, p, ι, Q)
    radiationforce = calculateradiationforce(mass, spin, massratio, initialvelnew, initialposnew, E, Lz, Eflux, Lflux, Qflux)
    inspirals = calculateinspiraleqns(mass, spin, radiationforce)

    tspan = (0, 12) # arbitrary value
    parameters = (inspirals,)
    prob = SecondOrderODEProblem(inspiraleqns!, initialvelnew, initialposnew, tspan, parameters)
    inspiralsoln = solve(prob, RK4(), abstol=1e-10, reltol=1e-10)

    # Save inspiral orbit
    tins = inspiralsoln[5, :]
    rins = inspiralsoln[6, :]
    θins = inspiralsoln[7, :]
    ϕins = inspiralsoln[8, :]
    utins = inspiralsoln[1, :]
    urins = inspiralsoln[2, :]
    uθins = inspiralsoln[3, :]
    uϕins = inspiralsoln[4, :]

    for i = eachindex(tins)
        if (round(rins[i], digits=3) > round(1.05*horizon, digits=3))
            push!(tinspiral, tins[i])
            push!(rinspiral, rins[i])
            push!(θinspiral, θins[i])
            push!(ϕinspiral, ϕins[i])
            push!(utinspiral, utins[i])
            push!(urinspiral, urins[i])
            push!(uθinspiral, uθins[i])
            push!(uϕinspiral, uϕins[i])

            # Calculate and save acceleration
            tτdotdot = inspirals[1](utins[i], urins[i], uθins[i], uϕins[i], 
                                    tins[i], rins[i], θins[i], ϕins[i])

            rτdotdot = inspirals[2](utins[i], urins[i], uθins[i], uϕins[i], 
                                    tins[i], rins[i], θins[i], ϕins[i])

            θτdotdot = inspirals[3](utins[i], urins[i], uθins[i], uϕins[i], 
                                    tins[i], rins[i], θins[i], ϕins[i])

            ϕτdotdot = inspirals[4](utins[i], urins[i], uθins[i], uϕins[i], 
                                    tins[i], rins[i], θins[i], ϕins[i])
                        
            push!(utdotinspiral, tτdotdot)
            push!(urdotinspiral, rτdotdot)
            push!(uθdotinspiral, θτdotdot)
            push!(uϕdotinspiral, ϕτdotdot)
        end
    end

    if saveepιfluxes
        push!(tfluxes, tins[length(inspiralsoln.t)])
        push!(einspiral, e)
        push!(pinspiral, p)
        push!(ιinspiral, ι)
        push!(Efluxes, Eflux)
        push!(Lfluxes, Lflux)
        push!(Qfluxes, Qflux)
    end

    prevsoln = inspiralsoln
    # previous final conditions -> new initial conditions
    global initialvelnew, initialposnew = calculatelastvelocityposition(prevsoln)
    global e, p, ι = calculateepι(prevsoln)
    global E, Lz, Q = calculateELQfromvelpos(mass, spin, littlemass, initialvelnew, initialposnew)
    # global E, Lz, Q = calculateELQfromepι(mass, spin, e, p, ι)
end

# Interpolate inspiral orbit
numpoints = length(tinspiral)

# Interpolation functions
interpolater = linear_interpolation(tinspiral, rinspiral)
interpolateθ = linear_interpolation(tinspiral, θinspiral)
interpolateϕ = linear_interpolation(tinspiral, ϕinspiral)
interpolateut = linear_interpolation(tinspiral, utinspiral)
interpolateur = linear_interpolation(tinspiral, urinspiral)
interpolateuθ = linear_interpolation(tinspiral, uθinspiral)
interpolateuϕ = linear_interpolation(tinspiral, uϕinspiral)
interpolateutdot = linear_interpolation(tinspiral, utdotinspiral)
interpolateurdot = linear_interpolation(tinspiral, urdotinspiral)
interpolateuθdot = linear_interpolation(tinspiral, uθdotinspiral)
interpolateuϕdot = linear_interpolation(tinspiral, uϕdotinspiral)

titp = range(tinspiral[1], tinspiral[numpoints], numpoints)
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

# Calculate waveform

# Parameters
distance = convertpctonaturalunits(mass, 5e9) # 5e9 pc
bigtheta = π/4
bigphi = 0

# To save waveform components
hPluses = Vector{Float64}()
hCrosses = Vector{Float64}()

x, y, z = calculatexyz(ritp, θitp, ϕitp)
xτdot, yτdot, zτdot = calculatexyzτdot(ritp, θitp, ϕitp, uritp, uθitp, uϕitp)
xτdotdot, yτdotdot, zτdotdot = calculatexyzτdotdot(ritp, θitp, ϕitp, 
                                                uritp, uθitp, uϕitp,
                                                urdotitp, uθdotitp, uϕdotitp)

xtdot, ytdot, ztdot = calculatexyztdot(xτdot, yτdot, zτdot, utinspiral)
xtdotdot, ytdotdot, ztdotdot = calculatexyztdotdot(xτdotdot, yτdotdot, zτdotdot, utinspiral)

for i = 1:length(titp)
    xyzarr = [x[i], y[i], z[i]]
    xyzdotarr = [xtdot[i], ytdot[i], ztdot[i]]
    xyzdotdotarr = [xtdotdot[i], ytdotdot[i], ztdotdot[i]]

    hbar = calculatehbar(massratio, distance, xyzarr, xyzdotarr, xyzdotdotarr)
    hPlus, hCross = calculatehcomponents(bigtheta, bigphi, hbar)

    push!(hPluses, hPlus)
    push!(hCrosses, hCross)
end

# Save inspiral data to CSV file
df1 = DataFrame(t = titp, r = ritp, θ = θitp, ϕ = ϕitp,
                tdot = utitp, rdot = uritp, θdot = uθitp, ϕdot = uϕitp,
                tdotdot = utdotitp, rdotdot = urdotitp, θdotdot = uθdotitp, ϕdotdot = uϕdotitp)
CSV.write("inspiral-orbit.csv", df1)

# Save orbital parameters
if saveepιfluxes
    df2 = DataFrame(t = tfluxes,
                    e = einspiral, p = pinspiral, ι = ιinspiral,
                    Eflux = Efluxes, Lflux = Lfluxes, Qflux = Qfluxes)
    CSV.write("inspiral-fluxes.csv", df2)
end

# Save waveform
df3 = DataFrame(plus = hPluses, cross = hCrosses)
CSV.write("waveform.csv", df3)