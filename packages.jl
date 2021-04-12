using Pkg

dependencies = [
    "PyCall",
    "ColorSchemes",
    "ColorTypes",
    "Images",
    "Luxor",
    "CSV",
    "DataFrames",
    "Optim",
    "JuMP",
    "Clp",
    "Cbc",
    "GLPK",
    "PackageCompiler",
    "MultivariateStats",
    "Glob",
    "ImageIO",
    "Ipopt"]

Pkg.add(dependencies)
