# puCCiNi.jl

# Installation

## Download Julia
Download Julia [https://julialang.org/downloads/](https://julialang.org/downloads/). 

Ensure that you can start Julia from the terminal

```bash
[user@host ~]$ julia
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.7.3 (2022-05-06)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |

julia>
```

## Download Files

Use git to clone repository

```bash
[user@host ~]$ git clone https://github.com/mdp-aerosol-group/puCCiNi.jl.git
Cloning into 'puCCiNi.jl'...
remote: Enumerating objects: 11, done.
remote: Counting objects: 100% (11/11), done.
remote: Compressing objects: 100% (9/9), done.
remote: Total 11 (delta 0), reused 7 (delta 0), pack-reused 0
Receiving objects: 100% (11/11), 15.72 KiB | 15.72 MiB/s, done.
```

## Download Dependencies

Navigate to directory and run ```julia --project```
```bash
[user@host ~]$ cd puCCiNi.jl/src
[user@host src/]$ julia --project
```

Run ```Pkg.instantiate()```. This may take some time...

```julia
julia> using Pkg	

julia> Pkg.instantiate()
  Installing known registries into `~/.julia`
    Updating registry at `~/.julia/registries/General.toml`
   Installed Calculus ───────────────────────── v0.5.1
   Installed JpegTurbo_jll ──────────────────── v2.1.2+0
   Installed x265_jll ───────────────────────── v3.5.0+0
   Installed TreeViews ──────────────────────── v0.3.0
   Installed SIMDDualNumbers ────────────────── v0.1.1
   Installed libfdk_aac_jll ─────────────────── v2.0.2+0
   Installed DifferentialEquations ──────────── v7.2.0
   Installed OffsetArrays ───────────────────── v1.12.7
   Installed Libmount_jll ───────────────────── v2.35.0+0
   Installed GR_jll ─────────────────────────── v0.66.0+0
   Installed HypergeometricFunctions ────────── v0.3.11
   Installed LERC_jll ───────────────────────── v3.0.0+1

OMITTED

  [29816b5a] + LibSSH2_jll
  [c8ffd9c3] + MbedTLS_jll
  [14a3606d] + MozillaCACerts_jll
  [4536629a] + OpenBLAS_jll
  [05823500] + OpenLibm_jll
  [bea87d4a] + SuiteSparse_jll
  [83775a58] + Zlib_jll
  [8e850b90] + libblastrampoline_jll
  [8e850ede] + nghttp2_jll
  [3f19e933] + p7zip_jll

julia> 
```

## Run Trajectories

```bash
[usr@host src/]$ time julia --project trajectories.jl

________________________________________________________
Executed in   46.71 secs    fish           external
   usr time   44.43 secs  234.00 micros   44.43 secs
   sys time    2.67 secs  109.00 micros    2.67 secs
```

The file produces ```.csv``` files as output. Each particle is given as a time series of ```z``` (height in mm), ```Dp``` (size in m) and ```v``` (terminal velocity in m/s).

Edit the input conditions on ```trajectories.jl``` to set up the simulation as desired.` 
