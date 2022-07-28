using Plots
using Plots.PlotMeasures
using DataFrames
using CSV

include("ChamberModel.jl")

T, p = 298.15, 1013e2           # Temperature [K] and Pressure [Pa]
t, n, zt = 20.1, 100, 0.01      # ODE integ. time [s], number of particles [#], chamber height [m]
Dg, σ, κ = 120e-9, 1.3, 0.6     # Mode diameter [m], geometric standard deviation [-], kappa
a, b = 0.002, zt * 0.9          # Range of chamber that particles are modeled [min, max]
seed = 704                      # Random seed                      
ss = [0.2, 0.5]                 # Array with supersaturations

mp = ChamberModel.ChamberProblem(T, p, t, n, zt, Dg, σ, a, b, κ, ss)
solutions = ChamberModel.integrateChamberProblem(mp, seed);

map(1:length(ss)) do i
    sol = solutions[i]

    t = 0.0:0.1:20.0            # Time step and time end for export
    zs = map(i -> 1000.0 * (hcat(sol[i](t)...))[2, :], 1:mp.n)
    Dps = map(i -> 2.0 * (hcat(sol[i](t)...))[1, :], 1:mp.n)
    v = map(Dp -> ChamberModel.vt.(Dp, mp.T, mp.p), Dps)

    df1 = DataFrame(zs, :auto)
    df2 = DataFrames.rename(df1, [i => Symbol("z$i") for i = 1:mp.n])
    df3 = DataFrame(Dps, :auto)
    df4 = DataFrames.rename(df3, [i => Symbol("D$i") for i = 1:mp.n])
    df5 = DataFrame(v, :auto)
    df6 = DataFrames.rename(df5, [i => Symbol("v$i") for i = 1:mp.n])

    str = "$(ss[i])"
    df = hcat(DataFrame(t = collect(t)), df2, df4, df6) |> CSV.write("out" * str * ".csv")
end
