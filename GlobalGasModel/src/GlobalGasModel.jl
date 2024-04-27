module GlobalGasModel

using JuMP
using XLSX
using CSV
using Parameters
using Dates
# using BilevelJuMP
using DataFrames
using Printf
using Base.Threads

export GGM_Parameters,
    GGM_Sets,
    GGM_Results,
    build_GGM,
    solve_GGM,
    get_GGM_inputs,
    get_target_values

const BCMA_TO_MCMD = 1000/365

function get_most_recent_file(dir, substring)
    scenariofiles =
        filter(filename -> occursin(substring, filename), readdir(dir))
    latest = argmax(mtime(joinpath.(dir, scenariofiles)))
    return joinpath(dir, scenariofiles[latest])
end

function set_optimizer_attributes(model, solver_options)
    for (option, value) in solver_options
        set_optimizer_attribute(model, option, value)
    end
end

function _print_info(step, status, time_elapsed; status_number = "", offset = 0)
    colwidth = Dict(
        :step => 35,
        :status => 50,
        :status_number => 15,
        :time_elapsed => 25,
    )
    if isa(status_number, String)
        @info string(
            lpad("", (offset)),
            rpad(step * "...", colwidth[:step]),
            rpad(status * " " * status_number, colwidth[:status]),
            rpad("", 5),
            rpad("Time elapsed in " * step * ":", colwidth[:step] + 15),
            rpad(
                Printf.@sprintf("%.2fs", time_elapsed),
                colwidth[:time_elapsed] + 15,
            ),
        )
    else
        @info string(
            lpad("", (offset)),
            rpad(step * "...", colwidth[:step]),
            rpad(
                status * lpad(
                    Printf.@sprintf("%.5e", status_number),
                    colwidth[:status_number],
                ),
                colwidth[:status],
            ),
            rpad("", 5),
            rpad("Time elapsed in " * step * ":", colwidth[:step] + 15),
            rpad(
                Printf.@sprintf("%.2fs", time_elapsed),
                colwidth[:time_elapsed] + 15,
            ),
        )
    end
end

@with_kw struct GGM_CalibrationTargets
    calib_lambda
    # Target_Q_P_country
    Target_Q_S_country
end

@with_kw mutable struct GGM_Parameters
    scen
    disc
    days_d
    cost_pl
    cost_pq
    cost_a
    cost_x
    inv_a
    inv_x
    inv_w
    int
    slp
    cour
    l_a
    l_i
    cap_a
    cap_p
    cap_x
    cap_w
    d_a_max
    d_x_max
    d_w_max
    sales_lim
end

@with_kw mutable struct GGM_Sets
    A
    Y
    N
    W
    R
    D
    N_c
    N_p
    T
    storage_idxs
    Q_P_idxs
    Q_S_idxs
    F_A_idxs
    F_I_idxs
    F_X_idxs
    D_A_idxs
    D_X_idxs
    D_W_idxs
    T_of_N_p
    T_of_N_c
    comb_T_and_N_p
    comb_T_and_N_c
    comb_T_and_N
    comb_T_and_A
    N_c_of_T_N
    A_s_of_T_N
    A_e_of_T_N
    W_of_N
    A_of_T
    T_of_A
    AV
    AP
    AL
    AR
end

@with_kw mutable struct GGM_Results
    result_D_A
    result_D_X
    result_D_W
    result_Q_P
    result_Q_S
    result_F_A
    result_F_I
    result_F_X
    result_P   
    result_Q_C
    result_Q_C_annual
    result_Q_S_annual
    result_Q_P_annual
    result_F_A_annual
end

include(joinpath(@__DIR__, "SubModules", "logging.jl"))
include(joinpath(@__DIR__, "SubModules", "data_load.jl"))
include(joinpath(@__DIR__, "SubModules", "Model.jl"))
include(joinpath(@__DIR__, "SubModules", "helper.jl"))

end # module
