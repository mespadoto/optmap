using PackageCompiler

if length(ARGS) == 1
    sysimage_path = ARGS[1]
else
    sysimage_path = "sys_optmap.so"
end

println("Saving image as $sysimage_path")

#Add to this list other solver packages as needed
packages_to_compile = [:MultivariateStats, :Optim, :JuMP, :Clp, :Cbc, :GLPK, :Ipopt, :PyCall, :ColorSchemes, :ColorTypes, :Images, :Luxor, :CSV, :DataFrames, :Glob]

create_sysimage(packages_to_compile,
                sysimage_path=sysimage_path, 
                precompile_execution_file="test_optmap.jl")
