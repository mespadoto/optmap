using JuMP
using Cbc

include("../OptMaps.jl")
using .OptMaps

function solve_model(solver, model_log, solver_options=Dict())
    model = Model(solver)

    for (k, v) in solver_options
        set_optimizer_attribute(model, k, v)
    end

    n_variables = 7
    n_constraints = 4

    @variable(model, 0 <= x[1:6] <= 10, Int)
    @variable(model, 0 <= y <= 1, Int)
    
    ocoeffs = [60 70 40 70 20 90]
    ccoeffs = [30 20 30 90 30 70]
    
    @constraint(model, sum((x[i]*ccoeffs[i] for i in 1:6)) <= 2000)
    @constraint(model, x[3] <= 10*x[4])
    @constraint(model, x[1] + x[2] >= 4*y)
    @constraint(model, x[5] + x[6] >= 4*(1 - y))

    constraints_f = []
    push!(constraints_f, (v) -> sum((v[i]*ccoeffs[i] for i in 1:6)) <= 2000)
    push!(constraints_f, (v) -> v[3] <= 10*v[4])
    push!(constraints_f, (v) -> v[1] + v[2] >= 4*v[7])
    push!(constraints_f, (v) -> v[5] + v[6] >= 4*(1 - v[7]))

    @objective(model, Max, sum(x[i]*ocoeffs[i] for i in 1:6))
    obj_f(v) = sum(v[i]*ocoeffs[i] for i in 1:6)

    starting_point = zeros(Float64, n_variables)
    solution_point = zeros(Float64, n_variables)

    ###############################################################
    #save model textual description
    open(model_log, "w") do f
        print(f, model)

        ###############################################################
        #solve model
        elapsed_time = @elapsed optimize!(model)

        println(f, "-----------------------------------------------")
        println(f, "Elapsed time: $elapsed_time")
        println(f, "-----------------------------------------------")
        println(f, "Status:    $(termination_status(model))")
        println(f, "Objective: $(objective_value(model))")
        println(f, "-----------------------------------------------")
        println(f, "Variables:")

        for i in 1:(n_variables-1)
            println(f, "$(x[i]): $(value(x[i]))")
            solution_point[i] = value(x[i])
        end

        println(f, "$(y): $(value(y))")
        solution_point[n_variables] = value(y)
    end

    variable_names = vcat(["x[$(string(i))]" for i in 1:(n_variables-1)], ["y"])
    constraint_names = ["c$(string(c))" for c in 1:n_constraints]

    obj_val = objective_value(model)

    return obj_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val
end

#The datatype used in the tuples is used to determine if the variable is "real" or integer
variable_domains = [(0, 12), (0, 12), (0, 12), (0, 12), (0, 12), (0, 12), (0, 1)]

optimizers = [Cbc.Optimizer]
opt_names = ["cbc"]

problems = Dict()

for (opt, opt_name) in zip(optimizers, opt_names)
    objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val = solve_model(opt, "knapsack_$(opt_name).log")
    problems[opt_name] = (objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val)
end

#All functions are equal, get first to initialize OptMap
objective_f, constraints_f, variable_names, constraint_names, _, _, obj_val = problems[opt_names[1]]
O = OptMap(variable_names, variable_domains, objective_f; constraint_names=constraint_names, constraints_f=constraints_f, verbose=true)
create_map!(O)

for opt_name in opt_names
    _, _, _, _, starting_point, solution_point, _ = problems[opt_name]

    path_to_solution = zeros(Float64, length(variable_domains), 2)
    path_to_solution[:,1] = starting_point
    path_to_solution[:,2] = solution_point

    draw_solution!(O, path_to_solution)
    save_all(O, "knapsack_$(opt_name)", ".")
end

