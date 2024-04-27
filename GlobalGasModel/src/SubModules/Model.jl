function build_GGM(
    Sets::GGM_Sets,
    Params::GGM_Parameters,
    Optimizer;
    solver_options = Dict(),
    store_log = "",
    calib_int = Dict(),
    offset = 0, 
    objective_scaling = 1,
    constraint_scaling = 1,
)
    starttime = time()
    
    @Status "Building Global Gas Model" "Progress" = "Starting." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log
    

    @unpack_GGM_Sets Sets
    @unpack_GGM_Parameters Params

    if isempty(calib_int)
        for k in keys(int)
            calib_int[k] = 1
        end
    end

    GGM = Model(Optimizer)

    @Status "Building Global Gas Model" "Progress" = "Creating variables." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log

    @variables GGM begin
        D_A[D_A_idxs] >= 0                        # Arc Capacity expansion (mcm / year)
        D_X[D_X_idxs] >= 0                      # Stor Extr Capacity expansion (mcm / year)
        D_W[D_W_idxs] >= 0                      # Stor WG Capacity expansion (mcm / year)

        Q_P[Q_P_idxs] >= 0                   # Quantity produced by resource (t is auxiliary) (mcm/yr)
        Q_S[Q_S_idxs] >= 0                   # Quantity sold (mcm/yr)

        F_A[F_A_idxs] >= 0                   # Trader Arc Flow (mcm/yr)
        F_I[F_I_idxs] >= 0                   # Trader Stor Injection (mcm/yr)
        F_X[F_X_idxs] >= 0                   # Trader Stor Extraction (mcm/yr)
    end

    #-------------------------------------------------------------------------------
    #---------------------    OBJECTIVE
    #-------------------------------------------------------------------------------

    @Status "Building Global Gas Model" "Progress" = "Creating objective." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log

    @expression(
        GGM,
        TC_QP,
        sum(
            disc[y] *
            days_d[d] *
            (
                cost_pl[n, r, y] *
                Q_P[(t, n, r, d, y)] +
                0.5 *
                cost_pq[n, r, y] *
                Q_P[(t, n, r, d, y)]^2
            ) for (t, n, r, d, y) in Q_P_idxs
        )
    )

    @expression(
        GGM,
        TC_FA,
        sum(
            disc[y] * days_d[d] * cost_a[a, y] * F_A[(t, a, d, y)] for
            (t, a, d, y) in F_A_idxs
        )
    )

    @expression(
        GGM,
        TC_FX,
        sum(
            disc[y] * days_d[d] * cost_x[n, w, y] * F_X[(t, n, w, d, y)] for
            (t, n, w, d, y) in F_X_idxs
        )
    )

    @expression(
        GGM,
        TC_inv,
        sum(disc[y] * inv_a[a, y] * D_A[(a, y)] for (a,y) in D_A_idxs ) +
        sum(
            disc[y] * (inv_x[n, w, y] * D_X[(n, w, y)]) for (n, w, y) in D_X_idxs
        ) +
        sum(
            disc[y] * (inv_w[n, w, y] * D_W[(n, w, y)]) for (n, w, y) in D_W_idxs
        )
    )

    @expression(GGM, TC, TC_QP + TC_FA + TC_FX + TC_inv)

    @expression(
        GGM,
        MPA,
        0.5 * sum(
            disc[y] *
            days_d[d] *
            slp[n_c, d, y] *
            cour[t, n_c, y] *
            (Q_S[(t, n_c, d, y)])^2 for (t, n_c, d, y) in Q_S_idxs
        )
    )

    @expression(
        GGM,
        Q_C[n_c ∈ N_c, d ∈ D, y ∈ Y],
        sum(Q_S[(t2, n_c, d, y)] for t2 in T_of_N_c[n_c])
    )

    @expression(
        GGM,
        REV,
        sum(
            disc[y] *
            days_d[d] *
            (calib_int[(n_c, d, y)] * int[n_c, d, y] - slp[n_c, d, y] * Q_C[n_c, d, y]) *
            Q_S[(t, n_c, d, y)] for (t, n_c, d, y) in Q_S_idxs
        )
    )

    @expression(
        GGM,
        CS,
        0.5 * sum(
            disc[y] * days_d[d] * slp[n_c, d, y] * (Q_C[n_c, d, y])^2 for
            n_c in N_c, d in D, y in Y
        )
    )

    @expression(GGM, P[n_c ∈ N_c, d ∈ D, y ∈ Y], calib_int[(n_c, d, y)] * int[n_c, d, y] - slp[n_c, d, y] * Q_C[n_c, d, y])

    @expression(GGM, Q_C_annual[n_c ∈ N_c, y ∈ Y], sum(days_d[d] * Q_C[n_c, d, y] for d ∈ D)/1000 )
    @expression(GGM, Q_S_annual[n_c ∈ N_c, t ∈ T_of_N_c[n_c], y ∈ Y], sum(days_d[d] * Q_S[(t, n_c, d, y)] for d ∈ D)/1000 )
    @expression(GGM, Q_P_annual[n_p ∈ N_p, t ∈ T_of_N_p[n_p], y ∈ Y], sum(days_d[d] * Q_P[(t, n_p, r, d, y)] for r ∈ R, d ∈ D)/1000 )
    @expression(GGM, F_A_annual[a ∈ A, t ∈ T_of_A[a], y in Y], sum( days_d[d] * F_A[(t, a, d, y)] for d in D )/1000)

    @objective( GGM, Max, objective_scaling * (REV + CS - TC - MPA) )

    @Status "Building Global Gas Model" "Progress" = "Creating constraints." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log

    @constraints GGM begin
        #-------------------------------------------------------------------------------
        #----------             SUPPLIER MASS BALANCES
        # -------------------------------------------------------------------------------
        eq_mass_bal[t in T, n in N, d in D, y in Y],
        constraint_scaling * (
        sum(Q_P[(t, n, r, d, y)] for r in R if (t, n) in comb_T_and_N_p) 
        +sum(F_A[(t, a, d, y)] * (1 - l_a[a]) for a in A_e_of_T_N[t, n]) 
        +sum(F_X[(t, n, w, d, y)] for w in W_of_N[n])
        ) ==
        constraint_scaling * (
        sum(Q_S[(t, n_c, d, y)] for n_c in N_c_of_T_N[t,n]) 
        +sum(F_A[(t, a, d, y)] for a in A_s_of_T_N[t, n]) 
        +sum(F_I[(t, n, w, d, y)] for w in W_of_N[n])
        )
        
        eq_stor_cycle[(t, n, w, d_, y) in F_X_idxs],
        constraint_scaling * (sum(days_d[d] * F_I[(t, n, w, d, y)] for d in D) * (1 - l_i[n, w])) ==
        constraint_scaling * (sum(days_d[d] * F_X[(t, n, w, d, y)] for d in D))
        #-------------------------------------------------------------------------------
        #-----------------       CAPACITY RESTRICTIONS
        #-------------------------------------------------------------------------------
        eq_cap_a[(a,y) in D_A_idxs, d in D],
        constraint_scaling * (sum(F_A[(t, a, d, y)] for t in T_of_A[a])) <=
        constraint_scaling * (cap_a[a, y] + sum(D_A[(a, y2)] for y2 in Y if y2 < y))
        
        eq_cap_p[(t, n, r, d, y) in Q_P_idxs],
        constraint_scaling * (Q_P[(t, n, r, d, y)]) <= constraint_scaling * (cap_p[n, r, y])
        
        eq_cap_x[(n, w, y) in D_X_idxs, d in D],
        constraint_scaling * (sum(F_X[(t, n, w, d, y)] for t in T_of_N_c[n])) <=
        constraint_scaling * (cap_x[n, w, y] + sum(D_X[(n, w, y2)] for y2 in Y if y2 < y))
        
        eq_cap_w[(n, w, y) in D_W_idxs],
        constraint_scaling * (sum(
            days_d[d] * F_X[(t, n, w, d, y)] for t in T_of_N_c[n], d in D
        )) <= constraint_scaling * (cap_w[n, w, y] + sum(D_W[(n, w, y2)] for y2 in Y if y2 < y))
        # ------------------------------------------------------------------------------
        # --------------------  INVESTMENT RESTRICTIONS
        # ------------------------------------------------------------------------------
        eq_lim_a[(a,y) in D_A_idxs], constraint_scaling * (D_A[(a, y)] )<= constraint_scaling * (d_a_max[a, y])
        eq_lim_x[(n, w, y) in D_X_idxs], constraint_scaling * (D_X[(n, w, y)]) <= constraint_scaling * (d_x_max[n, w, y])
        eq_lim_w[(n, w, y) in D_W_idxs], constraint_scaling * (D_W[(n, w, y)]) <= constraint_scaling * (d_w_max[n, w, y])

        # ------------------------------------------------------------------------------
        # --------------------  SALES RESTRICTIONS
        # ------------------------------------------------------------------------------
        eq_lim_sales[ (t, n_c, y) in keys(sales_lim)], sum(days_d[d] * Q_S[(t, n_c, d, y)] for d in D) <= sales_lim[t,n_c,y]
    end

    @Status "Building Global Gas Model" "Progress" = "Setting optimizer attributes." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log


    set_optimizer_attributes(GGM, solver_options)

    @Status "Building Global Gas Model" "Progress" = "Finished." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log

    return GGM
end


function solve_GGM(GGM; store_result = "", store_log = "", offset = 0, )

    starttime = time()

    @Status "Solving Global Gas Model" "Progress" = "Starting." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log

    optimize!(GGM)

    if termination_status(GGM) in [OPTIMAL, LOCALLY_SOLVED, ALMOST_LOCALLY_SOLVED] 
        @Status "Solving Global Gas Model" "Progress" = "Solve Successful." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log
    else 
        @GGMError "Model did not solve correctly" "Termination status" = termination_status(GGM) "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log
    end

    @Status "Solving Global Gas Model" "Progress" = "Getting results." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log
    
    res_idx_by_symb = Dict(
        :result_D_A => (:Arc, :Year),
        :result_D_X => (:Node, Symbol("Storage Type"), :Year),
        :result_D_W => (:Node, Symbol("Storage Type"), :Year),
        :result_Q_P => (:Trader, :Node, :Resource, :Season, :Year),
        :result_Q_S => (:Trader, :Node, :Season, :Year),
        :result_F_A => (:Trader, :Arc, :Season, :Year),
        :result_F_I => (:Trader, :Node, Symbol("Storage Type"), :Season, :Year),
        :result_F_X => (:Trader, :Node, Symbol("Storage Type"), :Season, :Year),
        :result_P => [:Node, :Season, :Year],
        :result_Q_C => [:Node, :Season, :Year],
        :result_Q_C_annual => [:Node, :Year],
        :result_Q_S_annual => [:Node, :Trader, :Year],
        :result_Q_P_annual => [:Node, :Trader, :Year],
        :result_F_A_annual => [:Arc, :Trader, :Year],
    )
    result_D_A = _to_df(round.(value.(GGM[:D_A]), digits=5), res_idx_by_symb[:result_D_A])
    println("result_D_A")
    result_D_X = _to_df(round.(value.(GGM[:D_X]), digits=5), res_idx_by_symb[:result_D_X])
    println("result_D_X")
    result_D_W = _to_df(round.(value.(GGM[:D_W]), digits=5), res_idx_by_symb[:result_D_W])
    println("result_D_W")
    result_Q_P = _to_df(round.(value.(GGM[:Q_P]), digits=5), res_idx_by_symb[:result_Q_P])
    println("result_Q_P")
    result_Q_S = _to_df(round.(value.(GGM[:Q_S]), digits=5), res_idx_by_symb[:result_Q_S])
    println("result_Q_S")
    result_F_A = _to_df(round.(value.(GGM[:F_A]), digits=5), res_idx_by_symb[:result_F_A])
    println("result_F_A")
    result_F_I = _to_df(round.(value.(GGM[:F_I]), digits=5), res_idx_by_symb[:result_F_I])
    println("result_F_I")
    result_F_X = _to_df(round.(value.(GGM[:F_X]), digits=5), res_idx_by_symb[:result_F_X])
    println("result_F_X")
    result_P   = _to_df(round.(value.(GGM[:P]), digits=5), res_idx_by_symb[:result_P])
    println("result_P")
    result_Q_C = _to_df(round.(value.(GGM[:Q_C]), digits=5), res_idx_by_symb[:result_Q_C])
    result_Q_C_annual = _to_df(round.(value.(GGM[:Q_C_annual ]), digits=5), res_idx_by_symb[:result_Q_C_annual])
    println("result_Q_C")
    result_Q_S_annual = _to_df(round.(value.(GGM[:Q_S_annual ]), digits=5), res_idx_by_symb[:result_Q_S_annual])
    println("result_Q_S_annual")
    result_Q_P_annual = _to_df(round.(value.(GGM[:Q_P_annual ]), digits=5), res_idx_by_symb[:result_Q_P_annual])
    println("result_Q_P_annual")
    result_F_A_annual = _to_df(round.(value.(GGM[:F_A_annual ]), digits=5), res_idx_by_symb[:result_F_A_annual])
    println("result_F_A_annual")

    Results = @pack_GGM_Results

    if length(store_result) > 0
        @Status "Solving Global Gas Model" "Progress" = "Storing results." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log
        if store_result[end-4:end] == ".xlsx"
            _to_excel(Results, store_result, offset , store_log )
        else
            _to_csv(Results, store_result, offset , store_log )
        end
    end

    @Status "Solving Global Gas Model" "Progress" = "Solution procedure finsished." "Time elapsed" = temporal(time() - starttime) offset = offset logfile = store_log

    return GGM, Results
end