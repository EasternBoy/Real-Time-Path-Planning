"Parallel version"
function ADMM!(cen::central, robo::Vector{robot}, H::Int64, Neb::Vector{Vector{Int64}}, 
                    Jo::Float64, nv::Matrix{Float64}, nϕ::Matrix{Float64}; MAX_ITER = 100, thres = 1e-2)

    N = cen.N

    ξ     = Vector{opt_cl}(undef,N)
    w     = Vector{opt_cl}(undef,N)
    λ     = Matrix{opt_cl}(undef,N,N)
    v     = Matrix{Float64}(undef,H,N)
    ϕ     = Matrix{Float64}(undef,H,N)

    # Initial position setup for all vehicles
    InitX = Matrix{Float64}(undef, H, N)
    InitY = Matrix{Float64}(undef, H, N)
    for i in 1:N
        InitX[:,i] = robo[i].p[1]*ones(H)
        InitY[:,i] = robo[i].p[2]*ones(H)
    end


    for i in 1:N
        w[i] = opt_cl(InitX, InitY)
        for j in 1:N
            λ[i,j] = opt_cl(zeros(H,N), zeros(H,N))
        end
    end


    Jn = zeros(N)

    t0 = time_ns()
    for n in 1:MAX_ITER
        #println("Iterative step $l")
        # Parallel manner in each vehicle
        ρ = 1.0

        for i in 1:N
            ξ[i] = Local_Update!(robo, cen, i, H, w, λ[i,:], Neb, ρ, Jo)
        end
 
        for i in 1:N
            w[i], v[:,i], ϕ[:,i], Jn[i] = Constraint_Update!(robo[i], cen, i, H, w, ξ, λ[:,i], Neb, ρ, nv[:,i], nϕ[:,i])
        end

        # Update dual variable
        for i in 1:N
            for j in Neb[i]
                λ[i,j].Nx = λ[i,j].Nx + ρ*(ξ[i].Nx - w[j].Nx)
                λ[i,j].Ny = λ[i,j].Ny + ρ*(ξ[i].Ny - w[j].Ny)
            end
        end

        ter = 0
        
        for i in 1:N
            ter = ter + norm(vec(ξ[i].Nx - w[i].Nx))
            ter = ter + norm(vec(ξ[i].Ny - w[i].Ny))
        end

        ter = round(ter, digits = 4)
        # println(ter)
        
        if ter <= thres
            break
        end
    end

    t1 = time_ns()
    dt = (t1-t0)/1e9
    dt = round(dt/N; digits = 3)
    println("Average time for each: $dt (s)")
    println("")

    return w, v, ϕ
end



function Local_Update!(robo::Vector{robot}, cen::central, ii::Int64, H::Int64, w::Vector{opt_cl}, λ::Vector{opt_cl}, Neb::Vector{Vector{Int64}}, ρ::Float64, Jo::Float64)
    # Find pairs of robots that need collision avoidance
    N   = cen.N

    # opti = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    opti = Model(Ipopt.Optimizer)
    set_silent(opti)

    xF = robo[ii].pF[1]
    yF = robo[ii].pF[2]

    # Variables
    Nx = @variable(opti, [1:H, 1:N]) # slack variable
    Ny = @variable(opti, [1:H, 1:N]) # slack variable
    
    J = sum((Nx[h,ii] - xF)^2 + (Ny[h,ii] - yF)^2 for h in 1:H)
    # @variable(robo.opti, Jn)
    for j in ii+1:N
        if norm(robo[ii].p[1:2] - robo[j].p[1:2]) <= 1 && norm(robo[j].p[1:2] - robo[j].pF) > 0.1
            J = 0
            break
        end
    end

    # Objective function
    # α  = 0.5
    # J = sum((Nx[h,ii] - xF)^2 + (Ny[h,ii] - yF)^2 for h in 1:H)

    for j in Neb[ii]
        J = J + ρ/2*dot(vec(Nx - w[j].Nx + λ[j].Nx/ρ), vec(Nx - w[j].Nx + λ[j].Nx/ρ)) + ρ/2*dot(vec(Ny - w[j].Ny + λ[j].Ny/ρ), vec(Ny - w[j].Ny + λ[j].Ny/ρ))
    end

    @objective(opti, MOI.MIN_SENSE, J)
    optimize!(opti)

    return opt_cl(JuMP.value.(Nx),  JuMP.value.(Ny))
end





function Constraint_Update!(robo::robot, cen::central, ii::Int64, H::Int64, w::Vector{opt_cl}, ξ::Vector{opt_cl}, 
                            λ::Vector{opt_cl}, Neb::Vector{Vector{Int64}}, ρ::Float64, nv::Vector{Float64}, nϕ::Vector{Float64})
    # Find pairs of robots that need collision avoidance
    N   = cen.N

    # robo.opti = Model(() -> Gurobi.Optimizer(GUROBI_ENV))
    robo.opti = Model(Ipopt.Optimizer)
    set_silent(robo.opti)


    xF = robo.pF[1]
    yF = robo.pF[2]
    τ  = robo.τ 

    # Variables
    ϕ  = @variable(robo.opti, [1:H])
    v  = @variable(robo.opti, [1:H])
    Nx = @variable(robo.opti, [1:H, 1:N]) # slack variable
    Ny = @variable(robo.opti, [1:H, 1:N]) # slack variable
    
    # @variable(robo.opti, Jn)

    # Objective function
    # J  = 0
    # Jw = sum(v[h]^2 for h in 1:H)
    Jw = 0
    # @constraint(robo.opti, J  <= Jo)

    for j in Neb[ii]
        Jw = Jw + ρ/2*dot(vec(ξ[j].Nx - Nx + λ[j].Nx/ρ), vec(ξ[j].Nx - Nx + λ[j].Nx/ρ)) + ρ/2*dot(vec(ξ[j].Ny - Ny + λ[j].Ny/ρ), vec(ξ[j].Ny - Ny + λ[j].Ny/ρ))
    end
    @objective(robo.opti, MOI.MIN_SENSE, Jw)





    # Constraints
    @constraints(robo.opti, begin robo.v_min .<= v .<= robo.v_max end)
    for h in 1:H
        if h == 1
            @constraint(robo.opti, -0.1 <= ϕ[h] - robo.p[3] <= 0.1)
        else
            @constraint(robo.opti, -0.1 <= ϕ[h] - ϕ[h-1]    <= 0.1)
        end
    end
    
    for i in 1:N
        for h in 1:H
            @constraint(robo.opti, robo.x_min + robo.Ra .<= Nx[h,i] .<= robo.x_max - robo.Ra)
            @constraint(robo.opti, robo.y_min + robo.Ra .<= Ny[h,i] .<= robo.y_max - robo.Ra)
        end
    end

    # Linearized Constraints
    for h in 1:H
        if h == 1
            ϕ0 =  robo.p[3]
            v0 =  robo.v
            @constraint(robo.opti, Nx[h,ii] - robo.p[1]  - τ*(v0*cos(ϕ0) + cos(ϕ0)*(v[h] - v0) - sin(ϕ0)*v0*(ϕ[h] - ϕ0))  == 0)
            @constraint(robo.opti, Ny[h,ii] - robo.p[2]  - τ*(v0*sin(ϕ0) + sin(ϕ0)*(v[h] - v0) + cos(ϕ0)*v0*(ϕ[h] - ϕ0))  == 0)
        else
            @constraint(robo.opti, Nx[h,ii] - Nx[h-1,ii] - τ*(nv[h]*cos(nϕ[h]) + cos(nϕ[h])*(v[h] - nv[h]) - sin(nϕ[h])*nv[h]*(ϕ[h] - nϕ[h])) == 0)
            @constraint(robo.opti, Ny[h,ii] - Ny[h-1,ii] - τ*(nv[h]*sin(nϕ[h]) + sin(nϕ[h])*(v[h] - nv[h]) + cos(nϕ[h])*nv[h]*(ϕ[h] - nϕ[h])) == 0)
        end
    end

    #Trusted region
    for h in 1:H
        @constraint(robo.opti, -0.5 <= v[h] - nv[h] <= 0.5)
        @constraint(robo.opti, -0.2 <= ϕ[h] - nϕ[h] <= 0.2)
    end
    
    ϵ = 0.1

    # Collision avoidance
    for j in 1:N
        if j != ii
            for h in 1:H
                dis = sqrt((w[ii].Nx[h,ii] - w[ii].Nx[1,j])^2 + (w[ii].Ny[h,ii] - w[ii].Ny[1,j])^2)
                ex =       (w[ii].Nx[h,ii] - w[ii].Nx[1,j])/dis
                ey =       (w[ii].Ny[h,ii] - w[ii].Ny[1,j])/dis
                @constraint(robo.opti, (Nx[h,ii] - Nx[h,j])*ex + (Ny[h,ii] - Ny[h,j])*ey >= 2*robo.Ra + ϵ)
            end
        end
    end


    # Now solve and get results
    status = JuMP.optimize!(robo.opti)

    Jn = sum((JuMP.value.(Nx)[h,ii] - xF)^2 + (JuMP.value.(Ny)[h,ii] - yF)^2 for h in 1:H)


    return opt_cl(JuMP.value.(Nx),  JuMP.value.(Ny)), JuMP.value.(v), JuMP.value.(ϕ), Jn
end
