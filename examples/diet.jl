using JuMP
using Clp
using Cbc
using GLPK

include("../OptMaps.jl")
using .OptMaps

function solve_model(solver, model_log, solver_options=Dict())
    model = Model(solver)

    for (k, v) in solver_options
        set_optimizer_attribute(model, k, v)
    end

    #sets
    foods = [:carrots, :potatoes, :bread, :cheese]
    resources = [:calories, :fat, :protein, :carbohydrate]

    #parameters
    cost_per_food = [0.14 0.4 0.3 0.75]
    rsc_per_food = [23.0  171.0  65.0  112.0;
                     0.1    0.2   0.0    9.3;
                     0.6    3.7   2.2    7.0;
                     6.0   30.0  13.0    5.0]

    demand_per_rsc = [2000 30 200 250]

    n_foods = length(foods)
    n_resources =length(resources)

    n_variables = n_foods
    n_constraints = n_resources

    ###############################################################
    #decision variables
    @variable(model, amt_per_food[1:n_foods] >= 0)

    starting_point = zeros(Float64, n_foods)
    solution_point = zeros(Float64, n_foods)

    for i in 1:n_foods
        v = start_value(amt_per_food[i])
        starting_point[i] = v === nothing ? 0.0 : v
    end

    ###############################################################
    #objective

    #regular function
    total_cost_f(amt_per_food) = sum([cost_per_food[f] * amt_per_food[f] for f in 1:n_foods])

    #jump expression
    total_cost = @expression(model, total_cost, sum(cost_per_food[f] * amt_per_food[f] for f in 1:n_foods))
    @objective(model, Min, total_cost)

    ###############################################################
    #constraints: fulfill min nutrition demand (demand_per_rsc)

    #regular functions
    constraints_f = []

    for r in 1:n_resources
        push!(constraints_f, (amt_per_food) -> sum([amt_per_food[f]*rsc_per_food[r,f]  for f in 1:n_foods]) >= demand_per_rsc[r])
    end

    #jump expressions
    @constraint(model, dmd[r = 1:n_resources], sum((amt_per_food[f]*rsc_per_food[r,f]  for f in 1:n_foods)) >= demand_per_rsc[r])

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

        for fd in 1:n_foods
            println(f, "$(foods[fd]): $(value(amt_per_food[fd]))")
            solution_point[fd] = value(amt_per_food[fd])
        end
        # println(f, "-----------------------------------------------")
        # println(f, "Shadow prices:")
        # for r in 1:n_resources
        #     println(f, "$r - $(resources[r]):  $(shadow_price(dmd[r]))")
        # end
    end

    variable_names = ["amt_per_food[$(String(f))]" for f in foods]
    constraint_names = ["dmd[$(String(r))]" for r in resources]

    obj_val = objective_value(model)

    return total_cost_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val
end

#The datatype used in the tuples is used to determine if the variable is "real" or integer
variable_domains = [(0.0,30.0), (0.0,30.0), (0.0,30.0), (0.0,30.0)]

optimizers = [Clp.Optimizer, GLPK.Optimizer, Cbc.Optimizer]
opt_names = ["clp", "glpk", "cbc"]

problems = Dict()

for (opt, opt_name) in zip(optimizers, opt_names)
    objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val = solve_model(opt, "diet_$(opt_name).log")
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
    save_all(O, "diet_$(opt_name)", ".")
end

