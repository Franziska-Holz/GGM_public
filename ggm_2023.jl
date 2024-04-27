using GlobalGasModel
using Dates
using Gurobi

# scen = "STEPS_NENO"
scen = ARGS[1]
Years = 2020:5:2060

logfile = joinpath(@__DIR__, "results", "logs", "warning_logs", "GGM_log_run_$(now())_$(scen).log")
data_file = joinpath(@__DIR__, "data_2023")
# transport_check_file = joinpath(@__DIR__, "results", "inputs", "GGM_transport_values_$(now())_$(scen).xlsx")
parameter_results_file = joinpath(@__DIR__, "results", "inputs", "GGM_parameter_values_$(now())_$(scen)")
results_file = joinpath(@__DIR__, "results", "outputs", "GGM_results_$(now())_$(scen)")

sets, params = get_GGM_inputs(
    scen,
    Years;
    data_file_path = data_file,
    # store_check_transport = transport_check_file,
    store_parameters = parameter_results_file,
    store_log = logfile,
);

GGM = build_GGM(
    sets,
    params,
    Gurobi.Optimizer;
    store_log = logfile,
    solver_options = Dict("NumericFocus" => 1)
    # objective_scaling = 1e-7,
    # constraint_scaling = 1e-7,
);

GGM, Results = solve_GGM(GGM;
    store_result = results_file,
    store_log = logfile,
);

