function _to_df(dict::Dict{String,V}, dim_name, value_name=:Value) where V
    ka = String[]
    va = V[]
    for (k,v) in dict
        push!(ka,k)
        push!(va,round(v,digits=5))
    end
    df = DataFrame([ka, va], [dim_name[1], value_name])
    return df
end

function _to_df(dict::Dict, dim_names, value_name=:Value)
    df = DataFrame([merge(NamedTuple{dim_names}(k), NamedTuple{(value_name,)}(round(v,digits=5))) for (k,v) in dict])
end

function _to_df(var::JuMP.Containers.DenseAxisArray, dim_names, value_name=:Value)
    if length(var.axes) == 1
        entries = var.axes[1]
        df = DataFrame([name => [entries[n][i] for n in eachindex(entries)] for (i,name) in enumerate(dim_names)])
        df[!,value_name] = var.data
    else
        table = Containers.rowtable(var; header = vcat(dim_names,value_name))
        df = DataFrames.DataFrame(table)
    end
    return df
end

function _to_df(var::JuMP.Containers.SparseAxisArray, dim_names, value_name=:Value)
    table = Containers.rowtable(var; header = vcat(dim_names,value_name))
    df = DataFrames.DataFrame(table)
    return df
end

function _to_excel(Params::GGM_Parameters, inputs_path, store_log; offset = 0)

    @unpack_GGM_Parameters Params
    par_by_symb = Dict(:disc => disc, :days_d => days_d, :cost_pl => cost_pl, :cost_pq => cost_pq, :cost_a => cost_a, :cost_x => cost_x, :inv_a => inv_a, :inv_x => inv_x, :inv_w => inv_w, :int => int, :slp => slp, :cour => cour, :l_a => l_a, :l_i => l_i, :cap_a => cap_a, :cap_p => cap_p, :cap_x => cap_x, :cap_w => cap_w, :d_a_max => d_a_max, :d_x_max => d_x_max, :d_w_max => d_w_max)
    par_idx_by_symb = Dict(:disc => (:Year,), :days_d => (:Season,), :cost_pl => (:Node,:Resource,:Year), :cost_pq => (:Node,:Resource,:Year), :cost_a => (:Arc, :Year), :cost_x => (:Node, Symbol("Storage Type"),:Year), :inv_a => (:Arc,:Year), :inv_x => (:Node,Symbol("Storage Type"),:Year), :inv_w => (:Node,Symbol("Storage Type"),:Year), :int => (:Node,:Season,:Year), :slp => (:Node,:Season,:Year), :cour => (:Trader,:Node,:Year), :l_a => (:Arc,), :l_i => (:Node,Symbol("Storage Type")), :cap_a => (:Arc,:Year), :cap_p => (:Node, :Resource,:Year), :cap_x => (:Node,Symbol("Storage Type"),:Year), :cap_w => (:Node,Symbol("Storage Type"),:Year), :d_a_max => (:Arc,:Year), :d_x_max => (:Node,Symbol("Storage Type"),:Year), :d_w_max => (:Node,Symbol("Storage Type"),:Year))
    XLSX.openxlsx(inputs_path, mode="w") do xf
        for symb_param in [:disc, :days_d, :cost_pl, :cost_pq, :cost_a, :cost_x, :inv_a, :inv_x, :inv_w, :int, :slp, :cour, :l_a, :l_i, :cap_a, :cap_p, :cap_x, :cap_w, :d_a_max, :d_x_max, :d_w_max]
            df = _to_df(par_by_symb[symb_param], par_idx_by_symb[symb_param])
            XLSX.addsheet!(xf,String(symb_param))
            if !isempty(par_by_symb[symb_param])
                XLSX.writetable!(xf[String(symb_param)], df)
            else
                @Warning "Parameter encountered that does not contain values, the respective entry is left empty." "Left out parameter" = "$(string(symb_param))" offset = offset logfile = store_log
            end
        end
    end
end

function _to_excel(Results::GGM_Results, results_file = "", offset = 0, store_log = "")
    
    @unpack_GGM_Results Results
    res_by_symb = Dict(
    :result_D_A => result_D_A,
    :result_D_X => result_D_X,
    :result_D_W => result_D_W,
    :result_Q_P => result_Q_P,
    :result_Q_S => result_Q_S,
    :result_F_A => result_F_A,
    :result_F_I => result_F_I,
    :result_F_X => result_F_X,
    :result_P => result_P,
    :result_Q_C => result_Q_C,
    :result_Q_C_annual => result_Q_C_annual,
    :result_Q_S_annual => result_Q_S_annual,
    :result_Q_P_annual => result_Q_P_annual,
    :result_F_A_annual => result_F_A_annual,
    )
    
    XLSX.openxlsx(results_file, mode="w") do xf
        for (symb_res, df) in res_by_symb
            XLSX.addsheet!(xf,String(symb_res))
            if !isempty(df)
                XLSX.writetable!(xf[String(symb_res)], df)
            else
                @Warning "Missing result value encountered, left out" "Missing result" =  "$(symb_res)" offset = offset logfile = store_log
            end
        end
    end
end

function _to_csv(Params::GGM_Parameters, inputs_path, store_log; offset = 0)

    @unpack_GGM_Parameters Params
    par_by_symb = Dict(:disc => disc, :days_d => days_d, :cost_pl => cost_pl, :cost_pq => cost_pq, :cost_a => cost_a, :cost_x => cost_x, :inv_a => inv_a, :inv_x => inv_x, :inv_w => inv_w, :int => int, :slp => slp, :cour => cour, :l_a => l_a, :l_i => l_i, :cap_a => cap_a, :cap_p => cap_p, :cap_x => cap_x, :cap_w => cap_w, :d_a_max => d_a_max, :d_x_max => d_x_max, :d_w_max => d_w_max)
    par_idx_by_symb = Dict(:disc => (:Year,), :days_d => (:Season,), :cost_pl => (:Node,:Resource,:Year), :cost_pq => (:Node,:Resource,:Year), :cost_a => (:Arc, :Year), :cost_x => (:Node, Symbol("Storage Type"),:Year), :inv_a => (:Arc,:Year), :inv_x => (:Node,Symbol("Storage Type"),:Year), :inv_w => (:Node,Symbol("Storage Type"),:Year), :int => (:Node,:Season,:Year), :slp => (:Node,:Season,:Year), :cour => (:Trader,:Node,:Year), :l_a => (:Arc,), :l_i => (:Node,Symbol("Storage Type")), :cap_a => (:Arc,:Year), :cap_p => (:Node, :Resource,:Year), :cap_x => (:Node,Symbol("Storage Type"),:Year), :cap_w => (:Node,Symbol("Storage Type"),:Year), :d_a_max => (:Arc,:Year), :d_x_max => (:Node,Symbol("Storage Type"),:Year), :d_w_max => (:Node,Symbol("Storage Type"),:Year))
    
    mkpath(inputs_path)

    for symb_param in [:disc, :days_d, :cost_pl, :cost_pq, :cost_a, :cost_x, :inv_a, :inv_x, :inv_w, :int, :slp, :cour, :l_a, :l_i, :cap_a, :cap_p, :cap_x, :cap_w, :d_a_max, :d_x_max, :d_w_max]
        df = _to_df(par_by_symb[symb_param], par_idx_by_symb[symb_param])
        if !isempty(par_by_symb[symb_param])
            CSV.write(joinpath(inputs_path, string(symb_param,".csv") ), df)
        else
            @Warning "Parameter encountered that does not contain values, the respective entry is left empty." "Left out parameter" = "$(string(symb_param))" offset = offset logfile = store_log
        end
    end
end

function _to_csv(Results::GGM_Results, results_file = "", offset = 0, store_log = "")
    
    @unpack_GGM_Results Results
    res_by_symb = Dict(
    :result_D_A => result_D_A,
    :result_D_X => result_D_X,
    :result_D_W => result_D_W,
    :result_Q_P => result_Q_P,
    :result_Q_S => result_Q_S,
    :result_F_A => result_F_A,
    :result_F_I => result_F_I,
    :result_F_X => result_F_X,
    :result_P => result_P,
    :result_Q_C => result_Q_C,
    :result_Q_C_annual => result_Q_C_annual,
    :result_Q_S_annual => result_Q_S_annual,
    :result_Q_P_annual => result_Q_P_annual,
    :result_F_A_annual => result_F_A_annual,
    )
    
    mkpath(results_file)
    for (symb_res, df) in res_by_symb
        if !isempty(df)
            CSV.write(joinpath(results_file, string(symb_res,".csv") ), df)
        else
            @Warning "Missing result value encountered, left out" "Missing result" =  "$(symb_res)" offset = offset logfile = store_log
        end
    end
end