using Optim
using Luxor
using JuMP
using Clp
using Cbc
using GLPK
using Ipopt
using Glob
using Base.Filesystem

include("OptMaps.jl")
using .OptMaps

function solve_model(solver, model_log, solver_options=Dict())
    model = Model(solver)

    for (k, v) in solver_options
        set_optimizer_attribute(model, k, v)
    end

    #sets
    foods = [:brownie, :icecream, :cola, :cheesecake]
    resources = [:calories, :chocolate, :sugar, :fat]

    #parameters
    cost_per_food = [50 20 30 80]
    rsc_per_food = [2000 200 150 500;
                    3 2 0 0;
                    2 2 4 4;
                    2 4 1 5]

    demand_per_rsc = [645 8 10 8]

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
variable_domains = [(0.0,4.0), (0.0,4.0), (0.0,4.0), (0.0,4.0)]

optimizers = [GLPK.Optimizer, Ipopt.Optimizer, Cbc.Optimizer]
opt_names = ["glpk", "ipopt", "cbc"]

for (opt, opt_name) in zip(optimizers, opt_names)
    objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val = solve_model(opt, "test_optmap_$(opt_name).log")
end

opt = Clp.Optimizer
opt_name = "clp"

objective_f, constraints_f, variable_names, constraint_names, starting_point, solution_point, obj_val = solve_model(opt, "test_optmap_$(opt_name).log")

path_to_solution = zeros(Float64, length(variable_domains), 2)
path_to_solution[:,1] = starting_point
path_to_solution[:,2] = solution_point

O = OptMap(variable_names, variable_domains, objective_f; constraint_names=constraint_names, constraints_f=constraints_f, verbose=true)
create_map!(O)

draw_solution!(O, path_to_solution)
save_all(O, "test_optmap_$(opt_name)", ".")

test_files = glob("test_optmap*.*", ".")

#Precompile Optim
f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
initial_x = [0.0, 0.0]

solvers = [NelderMead, LBFGS, GradientDescent, SimulatedAnnealing, ConjugateGradient, Newton]

for solver in solvers
    result = optimize(f, initial_x, solver(), Optim.Options(store_trace = true, extended_trace=true, iterations=50000))

    if solver == NelderMead
        tr = Optim.trace(result)
        points = [state.metadata["centroid"] for state in tr]
    else
        points = Optim.x_trace(result)
    end
end
