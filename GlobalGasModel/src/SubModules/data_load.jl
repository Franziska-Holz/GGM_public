# using XLSX, DataFrames, Parameters
# const BCMA_TO_MCMD = 1000/365
function get_GGM_inputs(
    scen,
    Y;
    data_file_path = pwd(),
    store_parameters = "",
    store_log = "",
    store_check_transport = "",
    Inf_val = 1e+4, 
    NaN_val = 0,
    verbose = true,
    check_transport_cutoff_rel = 0.02,
    check_transport_cutoff_abs = 1e+6,
)
    starttime = time()
    @Status "Data Load Global Gas Model" Progress = "Reading in Data." "Time elapsed" = temporal(time() - starttime) logfile = store_log
    
    scens = split(scen,"_")
    # scens = [scen,scen]
    other_assumptions_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Other Assumptions", infer_eltypes = true))
    arcs_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Arcs", infer_eltypes = true))
    general_market_power_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "General Market Power", infer_eltypes = true))
    specific_market_power_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Specific Market Power", infer_eltypes = true))
    storages_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Storages", infer_eltypes = true))
    vessel_distances_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Vessel Distances", infer_eltypes = true))
    nodes_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Nodes", infer_eltypes = true))
    seasons_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Seasons", infer_eltypes = true))
    resources_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Resources", infer_eltypes = true))
    sales_limitations = DataFrame(XLSX.readtable(joinpath(data_file_path,"general_data_$(scens[2]).xlsx"), "Sales Restrictions", infer_eltypes = true))
    projected_quantities = DataFrame(XLSX.readtable(joinpath(data_file_path,"projected_data_$(scens[1]).xlsx"), "Reference Quantities", infer_eltypes = true))
    sectoral_consumption_distribution = DataFrame(XLSX.readtable(joinpath(data_file_path,"projected_data_$(scens[1]).xlsx"), "Sector Consumption Distribution", infer_eltypes = true))
    seasonal_consumption_distribution = DataFrame(XLSX.readtable(joinpath(data_file_path,"projected_data_$(scens[1]).xlsx"), "Season Consumption Distribution", infer_eltypes = true))
    consumption_sectors = DataFrame(XLSX.readtable(joinpath(data_file_path,"projected_data_$(scens[1]).xlsx"), "Consumption Sectors", infer_eltypes = true))
    calibrated_global_losses = DataFrame(XLSX.readtable(joinpath(data_file_path,"calibrated_data_$(scens[1]).xlsx"), "Global Losses", infer_eltypes = true))
    calibrated_production_cost_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"calibrated_data_$(scens[1]).xlsx"), "Production Cost Calibration", infer_eltypes = true))
    calibrated_production_capacity_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"calibrated_data_$(scens[1]).xlsx"), "Production Capacity Calibration", infer_eltypes = true))
    calibrated_price_data = DataFrame(XLSX.readtable(joinpath(data_file_path,"calibrated_data_$(scens[1]).xlsx"), "Price Calibration", infer_eltypes = true))
    
    @Status "Data Load Global Gas Model" Progress = "Data Files Read." "Time elapsed" = temporal(time() - starttime) logfile = store_log

    vessel_distances_data = DataFrame(
        "Arc" => vec(["Ship_$(nₛ)_$(nₑ)" for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
        "Start" => vec([nₛ for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
        "End" => vec([nₑ for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
        "Type" => vec(["Vessel" for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
        "Shipping Distance (1000 Seamiles)" => vec([vessel_distances_data[findfirst(==(nₛ), vessel_distances_data[!,"Starting Node"]), nₑ] for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
        "Calibrated Usage Cost"	=> vec([1.0 for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
        "Calibrated Investment Cost" => vec([1.0 for nₑ in setdiff(names(vessel_distances_data),["Starting Node"]), nₛ in vessel_distances_data[:, "Starting Node"]]),
    )
    append!(arcs_data, vessel_distances_data; cols=:union)

    A = arcs_data[:, :Arc]
    AP = arcs_data[findall(==("Pipeline"), arcs_data.Type), :Arc]
    AL = arcs_data[findall(==("Liquefaction"), arcs_data.Type), :Arc]
    AR = arcs_data[findall(==("Regasification"), arcs_data.Type), :Arc]
    AV = arcs_data[findall(==("Vessel"), arcs_data.Type), :Arc]
    # Aₛ = Dict(arcs_data[i,:Arc] => arcs_data[i, :Start] for i in 1:nrow(arcs_data))
    # Aₑ = Dict(arcs_data[i,:Arc] => arcs_data[i, :End] for i in 1:nrow(arcs_data))
    
    T = unique(nodes_data[nodes_data.Production .== true, :Country])
    D = seasons_data[:, :Season]
    N = nodes_data[:, :Node]
    K = consumption_sectors[:, "Consumption Sector"]
    W = unique(storages_data[:, :Type])
    R = resources_data[:, :Resource]
    N_c = nodes_data[findall(nodes_data.Consumption), :Node]
    N_p = nodes_data[findall(nodes_data.Production), :Node]
    # N_w = nodes_data[findall(nodes_data.Storage), :Node]
    T_of_N_p =  # n_p_acc_t_fast
        Dict(n => nodes_data[findall(==(n), nodes_data.Node), :Country] for n in N_p)
    domestic_N_of_T = 
        Dict(t => nodes_data[findall(==(t), nodes_data.Country), :Node] for t in T)
    N_l_of_t = 
        Dict(t => nodes_data[(nodes_data.Country .== t) .&& nodes_data.Liquefaction, :Node] for t in T)
    domestic_N_c_of_T = 
        Dict( t => nodes_data[(nodes_data.Country .== t) .&& nodes_data.Consumption, :Node] for t in T)

    comb_T_and_N_p = Set( (nodes_data[i, :Country], nodes_data[i, :Node]) for i in findall(nodes_data.Production) ) # t_acc_np
    W_of_N = Dict(n => storages_data[findall(==(n), storages_data.Node), :Type] for n in N) # w_n_fast
    
    T_of_A, A_of_T, comb_T_and_A = get_arcs_of_trader(arcs_data, T, A, domestic_N_of_T)       # t_acc_a_fast
    
    comb_T_and_N = get_combinations_trader_active_nodes(T, N_l_of_t, arcs_data, nodes_data)
    # A_s_of_T_N = Dict( (t,n) => [a for a in arcs_data[ findall(==(n), arcs_data.Start ) ,:Arc] if (t,a) in comb_T_and_A ] for (t,n) in comb_T_and_N ) # a_s_fast
    # A_e_of_T_N = Dict( (t,n) => [a for a in arcs_data[ findall(==(n), arcs_data.End ) ,:Arc] if (t,a) in comb_T_and_A ] for (t,n) in comb_T_and_N ) # a_s_fast
    A_s_of_T_N = Dict( (t,n) => [a for a in arcs_data[ findall(==(n), arcs_data.Start ) ,:Arc] if (t,a) in comb_T_and_A ] for t in T, n in N) # a_s_fast
    A_e_of_T_N = Dict( (t,n) => [a for a in arcs_data[ findall(==(n), arcs_data.End ) ,:Arc] if (t,a) in comb_T_and_A ] for t in T, n in N) # a_s_fast

    disc = Dict(y => 1 / ((1 + other_assumptions_data[1, "Discount Rate"])^(y - Y[1])) for y in Y)
    days_d = Dict(d => other_assumptions_data[1, d] for d in D)
    days_y = sum(days_d[d] for d in D) 
    if days_y != 365 
        @Warning "Unexpected Number of Days encountered" "Total days in year" = sum(days_d[d] for d in D) logfile = store_log
    end
    cost_inflator = Dict(y => (1 + other_assumptions_data[1, "Cost Inflator"])^(y - Y[1]) for y in Y)
    price_inflator = Dict(y => (1 + other_assumptions_data[1, "Price Inflator"])^(y - Y[1]) for y in Y)
    
    ref_prod = Dict( (n,y) => projected_quantities[findfirst(==(n), projected_quantities.Node), "Reference Production $y"] * BCMA_TO_MCMD for n in N_p, y in Y )
    cap_p = Dict( (n, r, y) => ref_prod[n,y] * calibrated_production_capacity_data[findfirst(==(n), calibrated_production_capacity_data.Node), "$r"] * calibrated_production_capacity_data[findfirst(==(n), calibrated_production_capacity_data.Node), "$y"] for n in N_p, r in R, y in Y)
    sectoral_price_elasticity = Dict(k => consumption_sectors[findfirst(==(k),consumption_sectors[!,"Consumption Sector"]), "Price Elasticity"] for k in K)

    cost_pl = Dict((n, r, y) => calibrated_production_cost_data[findfirst(==(n), calibrated_production_cost_data.Node), "Base Cost (EUR/kcm)"] * calibrated_production_cost_data[findfirst(==(n), calibrated_production_cost_data.Node), "c($r)"] * calibrated_production_cost_data[findfirst(==(n), calibrated_production_cost_data.Node), "$y"] * cost_inflator[y] for n in N_p, r in R, y in Y)
    cost_pq = Dict((n, r, y) => calibrated_production_cost_data[findfirst(==(n), calibrated_production_cost_data.Node), "Base Cost (EUR/kcm)"] * calibrated_production_cost_data[findfirst(==(n), calibrated_production_cost_data.Node), "q($r)"] * calibrated_production_cost_data[findfirst(==(n), calibrated_production_cost_data.Node), "$y"] * cost_inflator[y] / cap_p[n,r,y] for n in N_p, r in R, y in Y)
    
    l_a = merge(Dict(a => other_assumptions_data[1, "Base Loss Liquefaction (Fraction Lost)"] for a in AL), 
        Dict(a => other_assumptions_data[1, "Base Loss Regasification (Fraction Lost)"] for a in AR), 
        Dict(a => other_assumptions_data[1, "Base Loss Pipeline (Fraction Lost per 1000km)"] * arcs_data[findfirst(==(a), arcs_data.Arc), "Length (1000 km)"] for a in AP), 
        Dict(a => other_assumptions_data[1, "Base Loss Ship (Fraction Lost per 1000 Sea Miles)"] * arcs_data[findfirst(==(a), arcs_data.Arc), "Shipping Distance (1000 Seamiles)"] for a in AV), 
    )

    cost_a = merge( 
        Dict( (a,y) => other_assumptions_data[1, "Base Fee Liquefaction (EUR/kcm)"] * (1 - l_a[a]) * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Usage Cost"] * cost_inflator[y] for a in AL, y in Y), 
        Dict( (a,y) => other_assumptions_data[1, "Base Fee Regasification (EUR/kcm)"] * (1 - l_a[a]) * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Usage Cost"] * cost_inflator[y] for a in AR, y in Y), 
        Dict( (a,y) => other_assumptions_data[1, "Base Fee Pipeline (EUR/kcm per 1000km)"]  * (1 - l_a[a]) * (arcs_data[findfirst(==(a), arcs_data.Arc), "Length (1000 km)"] + max(0, other_assumptions_data[1, "Multiplication Factor for Offshore Pipeline Investment"]-1) * arcs_data[findfirst(==(a), arcs_data.Arc), "Offshore Part  (1000 km)"] ) * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Usage Cost"] * cost_inflator[y] for a in AP, y in Y), 
        Dict( (a,y) => other_assumptions_data[1, "Base Fee Ship (EUR/kcm per 1000 Sea Miles)"] * (1 - l_a[a]) * arcs_data[findfirst(==(a), arcs_data.Arc), "Shipping Distance (1000 Seamiles)"] * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Usage Cost"] * cost_inflator[y] for a in AV, y in Y), 
    )

    inv_a = merge( 
        Dict( (a,y) => other_assumptions_data[1, "Base Investment Cost Liquefaction (EUR/kcm/d)"] * (1 - l_a[a]) * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Investment Cost"] * cost_inflator[y] / Y.step * days_y / days_y for a in AL, y in Y), 
        Dict( (a,y) => other_assumptions_data[1, "Base Investment Cost Regasification (EUR/kcm/d)"] * (1 - l_a[a]) * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Investment Cost"] * cost_inflator[y] / Y.step * days_y / days_y for a in AR, y in Y), 
        Dict( (a,y) => other_assumptions_data[1, "Base Investment Cost Pipeline (EUR/kcm/1000km)"]  * (1 - l_a[a]) * (arcs_data[findfirst(==(a), arcs_data.Arc), "Length (1000 km)"] + max(0, other_assumptions_data[1, "Multiplication Factor for Offshore Pipeline Investment"]-1) * arcs_data[findfirst(==(a), arcs_data.Arc), "Offshore Part  (1000 km)"] ) * arcs_data[findfirst(==(a), arcs_data.Arc), "Calibrated Investment Cost"] * cost_inflator[y] / Y.step * days_y / days_y for a in AP, y in Y), 
    )

    cap_a = Dict( (a,y) => arcs_data[findfirst(==(a), arcs_data.Arc), "$y Capacity (bcma)"] / (1 - l_a[a]) * BCMA_TO_MCMD for a in union(AL,AR,AP), y in Y)
    
    d_a_max = Dict( (a, Y[1]) => arcs_data[findfirst(==(a), arcs_data.Arc), "Maximum Expansion First Period"] / (1 - l_a[a]) * BCMA_TO_MCMD for a in union(AL,AR,AP) )
    if length(Y) > 1
        d_a_max = merge(d_a_max , Dict( (a, Y[2]) => arcs_data[findfirst(==(a), arcs_data.Arc), "Maximum Expansion Second Period"] / (1 - l_a[a])  * BCMA_TO_MCMD for a in union(AL,AR,AP) ))
        if length(Y) > 2    
            d_a_max = merge(d_a_max ,Dict( (a, y) => arcs_data[findfirst(==(a), arcs_data.Arc), "Maximum Expansion Other Periods"] / (1 - l_a[a])  * BCMA_TO_MCMD for a in union(AL,AR,AP), y in Y[3:end]))
        end
    end

    cost_x = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => storages_data[i, "OPEX (USD/kcm)"] * cost_inflator[y] for i in 1:nrow(storages_data), y in Y)
    inv_x = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => other_assumptions_data[1, "Base Investment Storage Extraction (EUR/kcm/d)"] * storages_data[i, "Extraction Calibration Factor"] * cost_inflator[y] / Y.step * days_y / days_y  for i in 1:nrow(storages_data), y in Y)
    cap_x = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => storages_data[i, "Extraction Capacity $y"] for i in 1:nrow(storages_data), y in Y)
    d_x_max = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => storages_data[i, "Extraction Capacity Expansion Limit"] for i in 1:nrow(storages_data), y in Y)

    inv_w = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => other_assumptions_data[1, "Base Investment Storage Working Gas (EUR/kcm)"] * storages_data[i, "Working Gas Calibration Factor"] * cost_inflator[y] / Y.step * days_y / days_y for i in 1:nrow(storages_data), y in Y)
    cap_w = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => storages_data[i, "Working Gas Capacity $y"] for i in 1:nrow(storages_data), y in Y)
    d_w_max = Dict( (storages_data[i,:Node], storages_data[i,:Type], y) => storages_data[i, "Working Gas Capacity Expansion Limit"] for i in 1:nrow(storages_data), y in Y)

    l_i = Dict( (storages_data[i,:Node], storages_data[i,:Type]) => storages_data[i, "Loss"] for i in 1:nrow(storages_data), y in Y )

    T_of_N_c, comb_T_and_N_c = get_combinations_trader_consumption_nodes(domestic_N_c_of_T, T, N_c, N_l_of_t, arcs_data)
    
    N_c_of_T_N = Dict( (t,n) => (((t,n) in comb_T_and_N_c ) ? [n] : [] ) for t in T, n in N )    # n_c_n_fast
    

    cour = Dict(
        (t, n_c, y) => 0.0 for t in T, n_c in N_c, y in Y
    )
    for t in general_market_power_data[:,:Trader]
        for n in N_c, y in Y
            if n ∉ domestic_N_c_of_T[t]
                cour[t, n, y] = max( general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Export"] * general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Minium Market Power Ratio"], general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Export"] * general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Moderation Factor"] ^ ((y - Y[1]) / Y.step) )
            else
                cour[t, n, y] = max( general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Domestic"] * general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Minium Market Power Ratio"], general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Domestic"] * general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Moderation Factor"] ^ ((y - Y[1]) / Y.step) )
            end
        end
    end
    for t in specific_market_power_data[:,:Trader]
        for n in names(specific_market_power_data, Not(:Trader)), y in Y
            cour[t, n, y] = max( specific_market_power_data[findfirst(==(t), specific_market_power_data.Trader), n] * general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Minium Market Power Ratio"], specific_market_power_data[findfirst(==(t), specific_market_power_data.Trader), n] * general_market_power_data[findfirst(==(t), general_market_power_data.Trader), "Moderation Factor"] ^ ((y - Y[1]) / Y.step) )
        end
    end

    global_reference_production = Dict( y => sum( projected_quantities[:,"Reference Production $y"] )  for y in Y)
    global_reference_consumption = Dict( y => sum( projected_quantities[:,"Reference Consumption $y"] )  for y in Y)

    sectoral_reference_consumption = Dict( (n,k,d,y) => (BCMA_TO_MCMD * projected_quantities[findfirst(==(n),projected_quantities.Node),"Reference Consumption $y"] * sectoral_consumption_distribution[findfirst(==(n), sectoral_consumption_distribution.Node), "$k $y"] * seasonal_consumption_distribution[findfirst(==(n), seasonal_consumption_distribution.Node), "$d"] * (1-calibrated_global_losses[findfirst(==(y), calibrated_global_losses.Year), "Global Loss"]) * global_reference_production[y] / global_reference_consumption[y]  )  for n in N_c, k in K, d in D, y in Y )
    reference_price = Dict( (n,d,y) => calibrated_price_data[findfirst(==(n), calibrated_price_data.Node), "Base Price (EUR/kcm)"] * calibrated_price_data[findfirst(==(n), calibrated_price_data.Node), "$d"] * calibrated_price_data[findfirst(==(n), calibrated_price_data.Node), "Price Calibration $y"] * price_inflator[y] for n in N_c, d in D, y in Y )

    seasonal_sector_intercept = Dict( (n,k,d,y) => reference_price[n,d,y] * (1-1/sectoral_price_elasticity[k]) for n in N_c, k in K, d in D, y in Y)
    seasonal_sector_slope = Dict( (n,k,d,y) => -reference_price[n,d,y] / (sectoral_price_elasticity[k] * sectoral_reference_consumption[n,k,d,y] ) for n in N_c, k in K, d in D, y in Y)

    ### 
    print(projected_quantities)
    ###
    slp = Dict( (n,d,y) => 1/sum(1/seasonal_sector_slope[n,k,d,y] for k in K) for n in N_c, d in D, y in Y)
    int = Dict( (n,d,y) => slp[n,d,y] * sum( seasonal_sector_intercept[n,k,d,y] / seasonal_sector_slope[n,k,d,y] for k in K) for n in N_c, d in D, y in Y )

    storage_idxs = Set( (storages_data[i,:Node], storages_data[i,:Type]) for i in 1:nrow(storages_data))
    Q_P_idxs = Set( (t, n, r, d, y) for (t,n) in comb_T_and_N_p, r in R, d in D, y in Y )
    Q_S_idxs = Set( (t, n, d, y) for (t,n) in comb_T_and_N_c, d in D, y in Y )
    F_A_idxs = Set( (t, a, d, y) for (t,a) in comb_T_and_A, d in D, y in Y)
    F_I_idxs = Set( (t, n, w, d, y) for t in T, (n, w) in storage_idxs, d in D, y in Y )
    F_X_idxs = Set( (t, n, w, d, y) for t in T, (n, w) in storage_idxs, d in D, y in Y )
    D_A_idxs = Set( (a,y) for a in union(AP,AL,AR), y in Y)
    D_X_idxs = Set( (n,w,y) for (n,w) in storage_idxs, y in Y)
    D_W_idxs = Set( (n,w,y) for (n,w) in storage_idxs, y in Y)

    sales_lim = Dict( (row[:Trader], row[:Node], row[:Year]) => row[:Limit] for row in eachrow(sales_limitations) )

    Sets = @pack_GGM_Sets
    Parameters = @pack_GGM_Parameters

    @Status "Data Load Global Gas Model" Progress = "Data Precalculated." "Time elapsed" = temporal(time() - starttime) logfile = store_log
    @Status "Data Load Global Gas Model" Progress = "Checking for Unexpected Values." "Time elapsed" = temporal(time() - starttime) logfile = store_log

    check_parameter_values(Sets, Parameters, nodes_data, Y, store_log, Inf_val = Inf_val, NaN_val = NaN_val)
    d_a_max_bcma = Dict(k => v/BCMA_TO_MCMD for (k,v) in d_a_max)
    check_transport_capacities(Y, projected_quantities, d_a_max_bcma, arcs_data, store_log; check_transport_cutoff_rel = check_transport_cutoff_rel, check_transport_cutoff_abs = check_transport_cutoff_abs)
    if length(store_check_transport) > 0
        transport_to_excel(Y, store_check_transport, projected_quantities, arcs_data, d_a_max_bcma)
    end
    if length(store_parameters) > 0
        @Status "Data Load Global Gas Model" Progress = "Storing Parameter Values." "Time elapsed" = temporal(time() - starttime) logfile = store_log
        if store_parameters[end-4:end] == ".xlsx"
            _to_excel(Parameters, store_parameters, store_log)
        else
            _to_csv(Parameters, store_parameters, store_log)
        end
    end
    
    @Status "Data Load Global Gas Model" Progress = "Data Processing Finished." "Time elapsed" = temporal(time() - starttime) logfile = store_log

    return Sets, Parameters
end

function transport_to_excel(Y, store_check_transport, projected_quantities, arcs_data, d_a_max_bcma)
    excel_df = DataFrame(Node = String[], Year = Int[], Production = Float64[], Consumption = Float64[], Liquefaction = Float64[], Pipeline_Export = Float64[], Regasification = Float64[], Pipeline_Import = Float64[])
    for i in 1:nrow(projected_quantities)
        for y in Y
            push!(excel_df, 
            (projected_quantities[i,"Node"], 
            y, 
            projected_quantities[i,"Reference Production $y"], 
            projected_quantities[i,"Reference Consumption $y"], 
            sum(arcs_data[j, "$y Capacity (bcma)"] + sum(d_a_max_bcma[arcs_data[j, "Arc"],y2] for y2 in Y if y2 < y ; init=0) for j in findall((arcs_data.Start .== projected_quantities[i,"Node"]) .& (arcs_data.Type .== "Liquefaction") ); init=0),
            sum(arcs_data[j, "$y Capacity (bcma)"] + sum(d_a_max_bcma[arcs_data[j, "Arc"],y2] for y2 in Y if y2 < y ; init=0) for j in findall((arcs_data.Start .== projected_quantities[i,"Node"]) .& (arcs_data.Type .== "Pipeline") ); init=0),
            sum(arcs_data[j, "$y Capacity (bcma)"] + sum(d_a_max_bcma[arcs_data[j, "Arc"],y2] for y2 in Y if y2 < y ; init=0) for j in findall((arcs_data.End .== projected_quantities[i,"Node"]) .& (arcs_data.Type .== "Regasification") ); init=0),
            sum(arcs_data[j, "$y Capacity (bcma)"] + sum(d_a_max_bcma[arcs_data[j, "Arc"],y2] for y2 in Y if y2 < y ; init=0) for j in findall((arcs_data.End .== projected_quantities[i,"Node"]) .& (arcs_data.Type .== "Pipeline") ); init=0)
            ))
        end
    end
    XLSX.writetable(store_check_transport, excel_df)
end

function check_transport_capacities(Y, projected_quantities, d_a_max_bcma, arcs_data, logfile = ""; check_transport_cutoff_rel = 0.02, check_transport_cutoff_abs = 1e+6)
    for i in 1:nrow(projected_quantities)
        for y in Y
            if (projected_quantities[i,"Reference Production $y"] > (1+check_transport_cutoff_rel)*projected_quantities[i,"Reference Consumption $y"]) || (projected_quantities[i,"Reference Production $y"] > projected_quantities[i,"Reference Consumption $y"] + check_transport_cutoff_abs)
                potential_export_capacity = sum(arcs_data[j, "$y Capacity (bcma)"] + sum(d_a_max_bcma[arcs_data[j, "Arc"],y2] for y2 in Y if y2 < y ; init = 0) for j in findall(==(projected_quantities[i,"Node"]), arcs_data.Start); init = 0)
                insufficient_export_capacity = round(projected_quantities[i,"Reference Production $y"] - projected_quantities[i,"Reference Consumption $y"] - potential_export_capacity;digits = 5)
                if insufficient_export_capacity > 0
                    @CheckInfo "Insufficient export capacity detected" "Year" = y "Node" = projected_quantities[i,"Node"] "Missing export capacity [BCMA]" = insufficient_export_capacity logfile = logfile
                end
            elseif (projected_quantities[i,"Reference Production $y"] < (1-check_transport_cutoff_rel)*projected_quantities[i,"Reference Consumption $y"]) || (projected_quantities[i,"Reference Production $y"] < projected_quantities[i,"Reference Consumption $y"] - check_transport_cutoff_abs)
                potential_import_capacity = sum(arcs_data[j, "$y Capacity (bcma)"] + sum(d_a_max_bcma[arcs_data[j, "Arc"],y2] for y2 in Y if y2 < y ; init=0) for j in findall(==(projected_quantities[i,"Node"]), arcs_data.End); init = 0)
                insufficient_import_capacity = round(projected_quantities[i,"Reference Consumption $y"] - projected_quantities[i,"Reference Production $y"] - potential_import_capacity;digits = 5)
                if insufficient_import_capacity > 0
                    @CheckInfo "Insufficient import capacity detected" "Year" = y "Node" = projected_quantities[i,"Node"] "Missing import capacity [BCMA]" = insufficient_import_capacity logfile = logfile
                end
            end
        end
    end
end

function get_arcs_of_trader(arcs_data, T, A, domestic_N_of_T)
    pipelines = arcs_data[findall(==("Pipeline"), arcs_data.Type), :Arc]
    regasification_arcs = arcs_data[findall(==("Regasification"), arcs_data.Type), :Arc]
    shipping_arcs = arcs_data[findall(==("Vessel"), arcs_data.Type), :Arc]
    A_of_T = Dict{String, Vector{String}}()
    comb_T_and_A = Set{Tuple{String, String}}()
    for t in T
        liquefaction_arcs = arcs_data[intersect(findall( n-> n in domestic_N_of_T[t], arcs_data.Start ),findall(==("Liquefaction"), arcs_data.Type)), :Arc]
        if !isempty(liquefaction_arcs)
            A_of_T[t] = union(liquefaction_arcs, pipelines, regasification_arcs, shipping_arcs)
            for arc in A_of_T[t]
                push!(comb_T_and_A, (t,arc))
            end
        else
            A_of_T[t] = pipelines
            for arc in A_of_T[t]
                push!(comb_T_and_A, (t,arc))
            end
        end
    end

    T_of_A = Dict( a => [t for (t,_a) in comb_T_and_A if a==_a ] for a in A )

    return T_of_A, A_of_T, comb_T_and_A
end

function get_combinations_trader_consumption_nodes(domestic_N_c_of_T::Dict{TT,Vector{TN}}, T, N_c, N_l_of_t, arcs_data) where {TT,TN}
    comb_T_and_N_c = Set{Tuple{TT, TN}}()
    T_of_N_c = Dict(n => TT[] for n in N_c)

    pipeline_ending_nodes = arcs_data[findall(==("Pipeline"), arcs_data.Type), :End]
    regasification_ending_nodes = arcs_data[findall(==("Regasification"), arcs_data.Type), :End]
    for (t,consumption_nodes) in domestic_N_c_of_T
        for consumption_node in consumption_nodes
            push!(comb_T_and_N_c, (t, consumption_node))
        end
    end
    for t in T
        for pipeline_ending_node in pipeline_ending_nodes
            if pipeline_ending_node in N_c
                push!(comb_T_and_N_c, (t, pipeline_ending_node))
            end
        end
        if !isempty(N_l_of_t[t])
            for regasification_ending_node in regasification_ending_nodes
                push!(comb_T_and_N_c, (t, regasification_ending_node))
            end
        end
    end

    for (t,n) in comb_T_and_N_c
        push!(T_of_N_c[n],t)
    end
    
    return T_of_N_c, comb_T_and_N_c
end

function get_combinations_trader_active_nodes(T, N_l_of_t, arcs_data, nodes_data)
    comb_T_and_N = Set( (nodes_data[i,:Country], nodes_data[i,:Node]) for i in 1:nrow(nodes_data) if nodes_data[i,:Production] )
    pipeline_ending_nodes = arcs_data[findall(==("Pipeline"), arcs_data.Type), :End]
    pipeline_starting_nodes = arcs_data[findall(==("Pipeline"), arcs_data.Type), :Start]
    regasification_nodes = arcs_data[arcs_data.Type .== "Regasification", :Start]
    for t in T
        for n in pipeline_ending_nodes
            push!(comb_T_and_N, (t,n) )
        end
        for n in pipeline_starting_nodes
            push!(comb_T_and_N, (t,n) )
        end
        liquefaction_nodes = N_l_of_t[t]
        if !isempty( liquefaction_nodes )
            for n in union(regasification_nodes, liquefaction_nodes)
                push!(comb_T_and_N, (t,n) )
            end
        end
    end
    return comb_T_and_N
end

function check_parameter_values(Sets, Parameters, nodes_data, Y, store_log; Inf_val = 1e+4, NaN_val = 0)
    @unpack_GGM_Sets Sets
    @unpack_GGM_Parameters Parameters
    var_by_symb = Dict(:disc => disc, :days_d => days_d, :cost_pl => cost_pl, :cost_pq => cost_pq, :cost_a => cost_a, :cost_x => cost_x, :inv_a => inv_a, :inv_x => inv_x, :inv_w => inv_w, :int => int, :slp => slp, :cour => cour, :l_a => l_a, :l_i => l_i, :cap_a => cap_a, :cap_p => cap_p, :cap_x => cap_x, :cap_w => cap_w, :d_a_max => d_a_max, :d_x_max => d_x_max, :d_w_max => d_w_max)
    for symb_param in [:disc, :days_d, :cost_pl, :cost_pq, :cost_a, :cost_x, :inv_a, :inv_x, :inv_w, :int, :slp, :cour, :l_a, :l_i, :cap_a, :cap_p, :cap_x, :cap_w, :d_a_max, :d_x_max, :d_w_max]
        param = var_by_symb[symb_param]
        for (k,v) in param
            if ismissing(v)
                @Warning "Unexpected parameter value encountered" "Parameter" = symb_param "Key" = k "Value encountered" = "Missing" "Set to" = 0 "logfile" = store_log
                param[k] = 0
            elseif v == Inf
                @Warning "Unexpected parameter value encountered" "Parameter" = symb_param "Key" = k "Value encountered" = "Inf" "Set to" = Inf_val "logfile" = store_log
                param[k] = Inf_val
            elseif isnan(v)
                @Warning "Unexpected parameter value encountered" "Parameter" = symb_param "Key" = k "Value encountered" = "NaN" "Set to" = NaN_val "logfile" = store_log
                param[k] = NaN_val
            end
        end
    end

    # All liquefiers in the nodes set should have exogenous capacity, or can endogenously be expanded, or both
    for n_l in nodes_data[nodes_data.Liquefaction, :Node]
        a = "LIQ_"*n_l
        try 
            liq_cap_sum = sum(cap_a[a, y] for y in Y) + sum(d_a_max[a, y] for y in Y)
            if liq_cap_sum <= 0
                @CheckInfo "Insufficient liquefaction capacity detected" "Node" = n_l "Capacity" = 0 "logfile" = store_log
            end
        catch e
            for y in Y
                if !haskey(cap_a, (a, y))
                    @CheckInfo "Missing entry detected" "Parameter" = "cap_a" "Node" = n_l "Arc" = a "Year" = y "logfile" = store_log
                end
                if !haskey(d_a_max, (a, y))
                    @CheckInfo "Missing entry detected" "Parameter" = "d_a_max" "Node" = n_l "Arc" = a "Year" = y "logfile" = store_log
                end
            end
        end
    end

    # All liquefiers in the nodes set should have positive operational costs, and positive investment costs
    for n_l in nodes_data[nodes_data.Liquefaction, :Node]
        a = "LIQ_"*n_l
        for y in Y
            if !haskey(cost_a, (a,y))
                @CheckInfo "Missing entry detected" "Parameter" = "cost_a" "Node" = n_l "Arc" = a "Year" = y "logfile" = store_log
            elseif cost_a[a,y] <= 0
                @CheckInfo "Unexpected value encountered" "Parameter" = "cost_a" "Node" = n_l "Arc" = a "Year" = y "value" = cost_a[a,y] "logfile" = store_log
            end
            if !haskey(inv_a, (a,y))
                @CheckInfo "Missing entry detected" "Parameter" = "inv_a" "Node" = n_l "Arc" = a "Year" = y "logfile" = store_log
            elseif inv_a[a,y] <= 0
                @CheckInfo "Unexpected value encountered" "Parameter" = "inv_a" "Node" = n_l "Arc" = a "Year" = y "value" = inv_a[a,y] "logfile" = store_log
            end
        end
    end

    # All liquefiers in the nodes set should have a loss rate in the range 5%-15%
    for n_l in nodes_data[nodes_data.Liquefaction, :Node]
        a = "LIQ_"*n_l
        if l_a[a] < 0.05 || l_a[a] > 0.15
            @CheckInfo "Unexpected value encountered" "Parameter" = "l_a" "Node" = n_l "Arc" = a "value" = l_a[a] "logfile" = store_log
        end
    end

    # All liquefiers in the nodes set should show up in the shipping distances matrix.
    liq_dif = setdiff(AL , "LIQ_" .* nodes_data[nodes_data.Liquefaction, :Node])
    if !isempty(liq_dif)
        for i in liq_dif
            @CheckInfo "Liquefier not present in both shipping nodes and shipping distance matrix" "Liquefier" = i "logfile" = store_log
        end
    end
    
    # All regasifiers in the nodes set should have exogenous capacity, or can endogenously be expanded, or both
    for n_r in nodes_data[nodes_data.Regasification, :Node]
        a = "REG_"*n_r
        try 
            reg_cap_sum = sum(cap_a[a, y] for y in Y) + sum(d_a_max[a, y] for y in Y)
            if reg_cap_sum <= 0 
                @CheckInfo "Insufficient regasification capacity detected" "Node" = n_r "Capacity" = 0 "logfile" = store_log
            end
        catch e
            for y in Y
                if !haskey(cap_a, (a, y))
                    @CheckInfo "Missing entry detected" "Parameter" = "cap_a" "Node" = n_r "Arc" = a "Year" = y "logfile" = store_log
                end
                if !haskey(d_a_max, (a, y))
                    @CheckInfo "Missing entry detected" "Parameter" = "d_a_max" "Node" = n_r "Arc" = a "Year" = y "logfile" = store_log
                end
            end
        end
    end

    # All regasifiers in the nodes set should have positive operational costs, and positive investment costs
    for n_r in nodes_data[nodes_data.Regasification, :Node]
        a = "REG_"*n_r
        for y in Y
            if !haskey(cost_a, (a,y))
                @CheckInfo "Missing entry detected" "Parameter" = "cost_a" "Node" = n_r "Arc" = a "Year" = y "logfile" = store_log
            elseif cost_a[a,y] <= 0
                @CheckInfo "Unexpected value encountered" "Parameter" = "cost_a" "Node" = n_r "Arc" = a "Year" = y "value" = cost_a[a,y] "logfile" = store_log
            end
            if !haskey(inv_a, (a,y))
                @CheckInfo "Missing entry detected" "Parameter" = "inv_a" "Node" = n_r "Arc" = a "Year" = y "logfile" = store_log
            elseif inv_a[a,y] <= 0
                @CheckInfo "Unexpected value encountered" "Parameter" = "inv_a" "Node" = n_r "Arc" = a "Year" = y "value" = inv_a[a,y] "logfile" = store_log
            end
        end
    end

    # All regasifiers in the nodes set should have a loss rate in the range 0.5%-3%
    for n_r in nodes_data[nodes_data.Regasification, :Node]
        a = "REG_"*n_r
        if l_a[a] < 0.005 || l_a[a] > 0.03
            @CheckInfo "Unexpected value encountered" "Parameter" = "l_a" "Node" = n_r "Arc" = a "value" = l_a[a] "logfile" = store_log
        end
    end

    # All regasifiers in the nodes set should show up in the shipping distances matrix
    reg_dif = setdiff(AR , "REG_" .* nodes_data[nodes_data.Regasification, :Node])
    if !isempty(reg_dif)
        for i in reg_dif
            @CheckInfo "Regasifier not present in both shipping nodes and shipping distance matrix" "Regasifier" = i "logfile" = store_log
        end
    end
end





