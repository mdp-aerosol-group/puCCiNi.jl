module ChamberModel

using Distributions
using Random
using Interpolations
using DifferentialEquations
using StatsBase

struct ChamberProblem
    T::Float64          # Temperature [K]     
    p::Float64          # Pressure [Pa]
    t::Float64          # ODE Integration time [s]
    n::Int              # Number of particles
    zt::Float64         # chamber height [m]
    Dg::Float64         # Mode diameter [m]
    σ::Float64          # geometric standard deviation [-]
    a::Float64          # bottom distance where particles appear
    b::Float64          # top where particles appear
    κ::Float64          # hygroscopicity of aerosol
    ss::Vector{Float64} # suppersaturations
end

const ρw = 997.1        # Water denity [kg m-3]
const lv = 2.5e6        # Latent heat of vaporization [J kg-1]
const R = 8.314         # Universal gas constant [J K-1 mol-1] = [N m K-1 mol-1]
const Mw = 18.05e-3     # Molecular weight of water [kg mol-1]
const Dv = 2.11e-5      # Diffusion constant of water [m2 s-1]
const K = 0.02587       # Heat conductivity [J s-1 m-1 K-1]
const g = 9.81          # Acceleration due to gravity [m s-2]
const Rd = 287.15       # Specific gas constant of dry air [J kg K]
const A = 2.1e-9        # Kelvin diameter [m] 
const ρp = 1000.0       # Droplet density [kg m-3]

# Saturation vapor pressure, T in [C]
es(T::Float64) = 610.94 * exp(17.625 * (T - 273.15) / ((T - 273.15) + 243.04))

# Mean free path of air
λ(T::Float64, p::Float64) = 6.6e-8 * (101315.0 / p) * (T / 293.15)

# viscosity of air
η(T::Float64) = 1.83245e-5 * exp(1.5 * log(T / 296.1)) * (406.55) / (T + 110.4)

# Slip correction factor
Cc(T::Float64, p::Float64, d::Float64) =
    1.0 + λ(T, p) / d * (2.34 + 1.05 * exp.(-0.39 * d / λ(T, p)))

# density of a gas
ρg(T::Float64, p::Float64) = p / (Rd * T)

# Reynolds number
Re(v::Float64, d::Float64, T::Float64, p::Float64) = ρg(T, p) * v * d / η(T)

# Drag coefficient
Cd(Re::Float64) = (Re < 0.1) ? 24.0 / Re : 24.0 / Re .* (1.0 + 0.15 * Re^0.687)

# Kohler theory/supersaturation over the drop
sᵈʳᵒᵖ(D::Float64, Dd::Float64, κ::Float64) =
    (D^3.0 - Dd^3.0) / (D^3.0 - Dd^3.0 * (1.0 - κ)) * exp(A / D) - 1.0

# Steady-state supersaturation profile in the CCN chamber vs. heigh
function sᶜʰᵃᵐᵇᵉʳ(z::Float64, ΔT::Float64, T::Float64, zt::Float64)
    e = es(T) + (es(T + ΔT) - es(T)) * z / zt
    s = (e / es(T + ΔT * z / zt) - 1)
    return (s < 0.0) ? 0.0 : s
end

# Growth parameter 
function G(T::Float64)
    FD = (ρw * R * T) / (es(T) * Dv * Mw)
    FH = (lv * ρw) / (K * T) * ((lv * Mw) / (T * R) - 1.0)
    return 1 / (FD + FH)
end

# Terminal velocity
function vt(d::Float64, T::Float64, p::Float64)
    vts1 =
        (d, v) -> (4.0 * ρp * d * g * Cc(T, p, d) / (3.0 * Cd(Re(v, d, T, p)) * ρg(T, p)))^0.5
    vts2 = (d, v) -> (vts1(d, v) - v)^2 < 1e-30 ? v : vts2(d, vts1(d, v))
    return vts2(d, 1e-5)
end

# Droplet trajectory ODEs
function parameterizedODE!(du, u, p, t)
    r, z = u
    G, Dd, ΔT, param = p
    du[1] = G * (sᶜʰᵃᵐᵇᵉʳ(z, ΔT, param.T, param.zt) - sᵈʳᵒᵖ(2.0 * r, Dd, param.κ)) / r
    du[2] = -vt(2.0 * r, param.T, param.p)
end

function integrateChamberProblem(p::ChamberProblem, seed::Int)
    smax(ΔT) = 100.0 * sᶜʰᵃᵐᵇᵉʳ.(0.0:0.0001:p.zt, ΔT, p.T, p.zt) |> maximum
    delT = 0.0:0.01:10.0
    mss = map(smax, delT)
    itp = interpolate((mss,), delT, Gridded(Linear()))
    Random.seed!(seed)
    
    function condition(u, t, integrator)
        out = ((u[1] < 1e-26) || (u[2] < 1e-26)) ? 0.0 : 1.0
    end
    
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)

    solutions = map(itp(p.ss)) do ΔT
        Dd, z0 = rand(LogNormal(log(p.Dg), log(p.σ)), p.n), rand(Uniform(p.a, p.b), p.n)
        Dd = map(D -> (D < 30e-9) ? 30e-9 : D, Dd)
        map(1:p.n) do i
            problem = ODEProblem(
                parameterizedODE!,
                [Dd[i]; z0[i]],
                (0.0, p.t),
                [G(p.T), Dd[i], ΔT, p],
            )
            solve(problem, alg_hints = [:stiff], abstol = 1e-11, reltol = 1e-11, callback = cb)
        end
    end

    return solutions
end

end
