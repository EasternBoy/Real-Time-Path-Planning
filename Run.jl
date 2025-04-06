push!(LOAD_PATH, ".")

import Pkg
Pkg.activate(".")
# Pkg.add("Ipopt")
# Pkg.build("HDF5")
# Pkg.precompile()


using LinearAlgebra, JuMP, Gurobi, Ipopt, Plots
using Dates

# if !(@isdefined(GUROBI_ENV))
#     const GUROBI_ENV = Gurobi.Env()
# end

# const GUROBI_ENV = Gurobi.Env()

include("Structure.jl")
include("ADMM.jl")
include("Initiate.jl")


# Initiate 
# Central unit: number of vehicles N, number of obstacle n_obs and sampling time ts
# Vehicle: maximum sampling times
# Obstacle 


st        = 0.1  # Sample time
H         = 10   # Horizontal prediction
sim_times = 60

cen, robo, obs = Initiate2!(st)

N = cen.N

w     = Vector{opt_cl}(undef, N)
traj  = zeros(3, sim_times + 1, N)

for i in 1:N
    traj[1, 1, i] = robo[i].p[1]
    traj[2, 1, i] = robo[i].p[2]
    traj[3, 1, i] = robo[i].p[3]
end    

InitX = Matrix{Float64}(undef, H, N)
InitY = Matrix{Float64}(undef, H, N)
for i in 1:N
    InitX[:,i] = robo[i].p[1]*ones(H)
    InitY[:,i] = robo[i].p[2]*ones(H)
    w[i]       = opt_cl(InitX, InitY)
end

J = 1e9
#--------------------------#
#Network define
Neb = Vector{Vector{Int64}}(undef,N)
for i in 1:N
    if i == 1
        Neb[i] = [N,   i, i+1]
    elseif i == N
        Neb[i] = [i-1, i, 1]
    else
        Neb[i] = [i-1, i, i+1]
    end
end
#--------------------------#
PredSteps  = zeros(2, H, sim_times + 1, N)

v  = zeros(H,N)
nv = zeros(H,N)
ϕ  = zeros(H,N)
nϕ = zeros(H,N)

for i in 1:N
    ϕ[:,i]  = robo[i].p[3]*ones(H)
    nϕ[:,i] = robo[i].p[3]*ones(H)
end


for k in 1:sim_times
    println("")

    if k == 1
        global nv = v
        global nϕ = ϕ
    else
        for i in 1:N
            global nv[1,i] = robo[i].v
            global nϕ[1,i] = robo[i].p[3]
            for h in 2:H
                opt = Model(Ipopt.Optimizer)
                set_silent(opt)
                @variable(opt,    robo[i].v_min    <= vopt <= robo[i].v_max)
                @variable(opt,    nϕ[h-1,i] - 0.1  <= ϕopt <= nϕ[h-1,i] + 0.1)
                @NLobjective(opt, MOI.MIN_SENSE, (w[i].Nx[h,i] - w[i].Nx[h-1,i] - vopt*cos(ϕopt))^2 + (w[i].Ny[h,i] - w[i].Ny[h-1,i] - vopt*sin(ϕopt))^2)
                JuMP.optimize!(opt)
                global nv[h,i] = JuMP.value.(vopt)
                global nϕ[h,i] = JuMP.value.(ϕopt)
            end
        end
    end

    global w, v, ϕ = ADMM!(cen, robo, H, Neb, J, nv, nϕ)
    println("Time step $k: v = $(v[1,1]), ϕ = $(ϕ[1,1])")

    
    for i in 1:N
        global PredSteps[1,:,k+1,i] =  w[i].Nx[:,i]
        global PredSteps[2,:,k+1,i] =  w[i].Ny[:,i]
    

        robo[i].p[1] = robo[i].p[1] + st*v[1,i]*cos(ϕ[1,i])
        robo[i].p[2] = robo[i].p[2] + st*v[1,i]*sin(ϕ[1,i])
        robo[i].p[3] = ϕ[1,i]
        robo[i].v    = v[1,i]


        global traj[1, k+1, i]      =  robo[i].p[1]
        global traj[2, k+1, i]      =  robo[i].p[2]
        global traj[3, k+1, i]      =  robo[i].p[3]
    end
end

matwrite("traj.mat", Dict("traj" => traj))
matwrite("PredSteps.mat", Dict("PredSteps" => PredSteps))



function circleShape(x,y,r)
    θ = LinRange(0, 2*π, 100)
    return x .+ r*cos.(θ), y .+ r*sin.(θ)
end


color = ["blue", "green", "pink", "yellow"]

for k = 1:2:sim_times
    p = plot()
    for i in 1:cen.N
        plot!(circleShape(traj.[1,k,i], traj.[2,k,i], robo[i].Ra), 
                            seriestype = [:shape,], legend = false, fillcolor = color[i], fillalpha = 1.0, tickfont = "Arial")
    end
    # xlims!(0, 10)
    # ylims!(0, 10)
    display(p)
end
