using JuMP, Gurobi
import OSQP

function Initiate2!(ts::Float64)

    n_veh = 4
    n_obs = 3

    init_pos_robo = [ 4.   3.
                    ; 3.   4.
                    ; 3.   3.
                    ; 4.   4.
                    ]
    int_ang_robo = [0, 3π/2, π/2, π]

    dest = [  7.  7.
            ; 9.  9.
            ; 9.  7.
            ; 7.  9.
            ]

    T = 20*ones(Int64,n_veh)


    cen  = central(n_veh, n_obs, ts, init_pos_robo)
    robo = Vector{robot}(undef, n_veh)
    obs  = Vector{obstacle}(undef, n_obs)


    v_min =-1.0
    v_max = 1.0
    x_min = 0.5
    x_max = 15.0
    y_min = 0.5
    y_max = 15.0
    Rd    = 2.0
    Ra    = 0.3


    # matwrite("init_pos.mat", Dict("init_pos_robo" => init_pos_robo))


    solver = Gurobi.Optimizer

    for i in 1:n_veh
        robo[i] = robot(T[i], ts, v_min, v_max, x_min, x_max, y_min, y_max, Rd, Ra, [init_pos_robo[i,:]; int_ang_robo[i]], dest[i,:])
    end



    obs[1]  = obstacle(1, 4.,   7.,   5.,  1., π/2) #Rectangular
    obs[2]  = obstacle(1, 10.0, 11.,  3.,  1., 0.)  #Rectangular
    obs[3]  = obstacle(1, 11.0, 4.,   4.,  1., 0.)  #Rectangular

    return cen, robo, obs

end



function Initiate3!(ts::Float64)

    n_veh = 4
    n_obs = 3

    init_pos_robo = [ 3.   3.
                    ; 3.   8.
                    ; 8.   3.
                    ; 8.   8.
                    ]
    int_ang_robo = [0, 0, 0, 0]

    dest = [  8.  8.
            ; 8.  3.
            ; 3.  8.
            ; 3.  3.
            ]

    T = 20*ones(Int64,n_veh)


    cen  = central(n_veh, n_obs, ts, init_pos_robo)
    robo = Vector{robot}(undef, n_veh)
    obs  = Vector{obstacle}(undef, n_obs)


    v_min =-1.0
    v_max = 1.0
    x_min = 0.5
    x_max = 15.0
    y_min = 0.5
    y_max = 15.0
    Rd    = 2.0
    Ra    = 0.3


    matwrite("init_pos.mat", Dict("init_pos_robo" => init_pos_robo))


    solver = Gurobi.Optimizer

    for i in 1:n_veh
        robo[i] = robot(T[i], ts, v_min, v_max, x_min, x_max, y_min, y_max, Rd, Ra, [init_pos_robo[i,:]; int_ang_robo[i]], dest[i,:])
    end



    obs[1]  = obstacle(1, 4.,   7.,   5.,  1., π/2) #Rectangular
    obs[2]  = obstacle(1, 10.0, 11.,  3.,  1., 0.)  #Rectangular
    obs[3]  = obstacle(1, 11.0, 4.,   4.,  1., 0.)  #Rectangular

    return cen, robo, obs

end