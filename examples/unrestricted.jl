using ColorSchemes
using ColorTypes
using Images
using Optim
using Random

include("functions.jl")
include("../OptMaps.jl")
using .OptMaps

dimensions = [3, 5, 7, 10, 20]
all_funcs = Dict()

func_types = ["sphere", "rastrigin", "rosenbrock", "styblinskitang", "linear"]
funcs = [get_sphere_funcs, get_rastrigin_funcs, get_rosenbrock_funcs, get_styblinskitang_funcs, get_linear_funcs]

for (ftype, fs) in zip(func_types, funcs)
    all_funcs[ftype] = fs(dimensions)
end

solvers = [NelderMead, LBFGS, GradientDescent, SimulatedAnnealing, ConjugateGradient, Newton]
solver_names = ["NelderMead", "LBFGS", "GradientDescent", "SimulatedAnnealing", "ConjugateGradient", "Newton"]

starting_points = ["half", "ones", "two_and_half", "threes", "rand"]
starting_points_elt = [0.5, 1.0, 2.5, 3.0, 0.0]

Random.seed!(787)

for ftype in func_types
    fs = all_funcs[ftype]

    for dim in dimensions
        println("Function: $ftype, dimensions: $dim")

        variable_domains = convert(Array{Tuple{Float64,Float64},1}, [])
        variable_names = convert(Array{String,1}, [])

        for i in 1:dim
            push!(variable_domains, (-5.0,5.0))
            push!(variable_names, string(i))
        end

        O = OptMap(variable_names, variable_domains, fs[dim]; verbose=true)
        create_map!(O)

        for (starting_point, starting_point_elt) in zip(starting_points, starting_points_elt)
            if starting_point == "rand"
                initial_x = (rand(Float64, dim) .* 10.0) .- 5.0
            else
                initial_x = zeros(Float64, dim)
                initial_x = initial_x .+ starting_point_elt
            end

            for (solver, solver_name) in zip(solvers, solver_names)
                result = optimize(fs[dim], initial_x, solver(), Optim.Options(store_trace = true, extended_trace=true); autodiff = :forward)
                result_log = "$(ftype)_$(dim)_$(solver_name)_$(starting_point)_optmap.log"

                open(result_log, "w") do f
                    println(f, "$(ftype), $(dim) dims, $solver_name - starting from $starting_point")
                    println(f, "$(initial_x)")
                    print(f, result)
                end            

                println("#######################################################")
                println("$solver_name - $starting_point")

                if solver_name == "NelderMead"
                    tr = Optim.trace(result)
                    points = [state.metadata["centroid"] for state in tr]
                else
                    points = Optim.x_trace(result)
                end

                solution_path = zeros(Float64, length(variable_names), length(points))

                for i in 1:length(points)
                    for j in 1:length(variable_names)
                        solution_path[j, i] = points[i][j]
                    end
                end

                draw_solution!(O, solution_path)
                save_all(O, "$(ftype)_$(dim)_$(solver_name)_$(starting_point)", ".")
            end
        end

        O = nothing

    end
end
