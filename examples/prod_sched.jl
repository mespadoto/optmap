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
    xs = [1, 2, 3, 4]
    ys = [:skilled, :unskilled]
    resources = [:mach1, :mach2]

    #parameters
    unit_prices = [300 260 220 180]
    labor_cost_hour = [8 6]

    rsc_per_prod = [11 7 6 5;
                    4 6 5 4]

    hours_per_prod = [8 5 5 6;
                      7 8 7 4]

    resource_limits = [700 500]
    hour_limits = [600 650]

    n_products = length(xs)
    n_labor_types = length(ys)
    n_resources = length(resources)

    n_variables = n_products + n_labor_types
    n_constraints = n_resources

    ###############################################################
    #decision variables
    @variable(model, x[1:n_products] >= 0)
    @variable(model, y[1:n_labor_types] >= 0)

    starting_point = zeros(Float64, n_variables)
    solution_point = zeros(Float64, n_variables)

    ###############################################################
    #objective

    #regular function
    total_profit_f(v) = sum([v[i] * unit_prices[i] for i in 1:n_products]) - sum(v[l+4] * labor_cost_hour[l] for l in 1:n_labor_types)
    
    #jump expression
    total_profit = @expression(model, total_profit, sum(x[p] * unit_prices[p] for p in 1:n_products) - sum(y[l] * labor_cost_hour[l] for l in 1:n_labor_types))

    @objective(model, Max, total_profit)

    ###############################################################
    #constraints:

    #regular functions
    constraints_f = []

    for r in 1:n_resources
        push!(constraints_f, (v) -> sum([v[p] * rsc_per_prod[r,p] for p in 1:n_products]) <= resource_limits[r])
    end

    for l in 1:n_labor_types
        push!(constraints_f, (v) -> sum([v[p] * hours_per_prod[l,p] for p in 1:n_products]) <= v[l+4])
    end

    for l in 1:n_labor_types
        push!(constraints_f, (v) -> v[l+4] <= hour_limits[l])
    end

    #jump expressions
    @constraint(model, rsc[r = 1:n_resources], sum((x[p]*rsc_per_prod[r,p]  for p in 1:n_products)) <= resource_limits[r])
    @constraint(model, hours[l = 1:n_labor_types], sum((x[p]*hours_per_prod[l,p]  for p in 1:n_products)) <= y[l])
    @constraint(model, hours_lim[l = 1:n_labor_types], y[l] <= hour_limits[l])

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

        for p in 1:n_products
            println(f, "$(xs[p]): $(value(x[p]))")
            solution_point[p] = value(x[p])
        end

        for l in 1:n_labor_types
            println(f, "$(ys[l]): $(value(y[l]))")
            solution_point[l+4] = value(y[l])
        end

        # println(f, "-----------------------------------------------")
        # println(f, "Shadow prices:")
        # for r in 1:n_resources
        #     println(f, "$r - $(resources[r]):  $(shadow_price(dmd[r]))")
        # end
    end

    variable_names = vcat(["x[$(string(p))]" for p in xs], ["y[$(string(l))]" for l in ys])
    constraint_names = ["mach1", "mach2", "skilled", "unskilled", "skilled_tot", "unskilled_tot"]

    obj_val = objective_value(model)

    return total_profit_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val
end

#The datatype used in the tuples is used to determine if the variable is "real" or integer
variable_domains = [(-100.0,100.0), (-100.0,100.0), (-100.0,100.0), (-100.0,100.0), (0.0,1000.0), (0.0,1000.0)]

optimizers = [Clp.Optimizer, GLPK.Optimizer, Cbc.Optimizer]
opt_names = ["clp", "glpk", "cbc"]

problems = Dict()

for (opt, opt_name) in zip(optimizers, opt_names)
    objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val = solve_model(opt, "prod_sched_$(opt_name).log")
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
    save_all(O, "prod_sched_$(opt_name)", ".")
end

