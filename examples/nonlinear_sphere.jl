using JuMP
using Ipopt

include("../OptMaps.jl")
using .OptMaps

function solve_model(solver, model_log, solver_options=Dict())
    model = Model(solver)

    for (k, v) in solver_options
        set_optimizer_attribute(model, k, v)
    end

    n_variables = 10
    n_constraints = 1

    @variable(model, -3.0 <= x[1:n_variables] <= 5.0)
    @NLconstraint(model, sum((x[i]^2 for i in 1:n_variables)) >= 5.0)

    constraints_f = []
    push!(constraints_f, (v) -> sum(v.^2) >= 5.0)

    @NLobjective(model, Min, sum(x[i]^2 for i in 1:n_variables))
    obj_f(v) = sum(v.^2)

    starting_point = zeros(Float64, n_variables) .+ 3.0
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

        for i in 1:n_variables
            println(f, "$(x[i]): $(value(x[i]))")
            solution_point[i] = value(x[i])
        end
    end

    variable_names = vcat(["x[$(string(i))]" for i in 1:n_variables])
    constraint_names = ["c1"]

    obj_val = objective_value(model)

    return obj_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val
end

optimizers = [Ipopt.Optimizer]
opt_names = ["ipopt"]

problems = Dict()

for (opt, opt_name) in zip(optimizers, opt_names)
    objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val = solve_model(opt, "sphere_$(opt_name).log")
    problems[opt_name] = (objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val)
end

#All functions are equal, get first to initialize OptMap
objective_f, constraints_f, variable_names, constraint_names, _, _, obj_val = problems[opt_names[1]]
variable_domains = convert(Array{Tuple{Float64,Float64},1}, [])

for i in 1:length(variable_names)
    push!(variable_domains, (-4.0,5.0))
end

O = OptMap(variable_names, variable_domains, objective_f; constraint_names=constraint_names, constraints_f=constraints_f, verbose=true)
create_map!(O)

for opt_name in opt_names
    _, _, _, _, starting_point, solution_point, _ = problems[opt_name]

    path_to_solution = zeros(Float64, length(variable_domains), 2)
    path_to_solution[:,1] = starting_point
    path_to_solution[:,2] = solution_point

    draw_solution!(O, path_to_solution; line_width=3.0)
    save_all(O, "sphere_$(opt_name)", ".")
end

