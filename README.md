# OptMap

This is the repository containing the source code for the paper titled "OptMap: Using Dense Maps for Visualizing Multidimensional Optimization Problems", by M. Espadoto, F. C. M. Rodrigues, N. S. T. Hirata and A. C. Telea, presented at IVAPP 2021.

## Julia Setup

- Install Julia 1.5.x from [here](https://julialang.org/downloads/)
- Install the required packages:

```
$ julia packages.jl
```

- Create a precompiled sysimage with the required packages:

```
$ julia precomp.jl
```

## Run Examples

- Check out the example scripts in the `examples` folder to see how OptMap is used. To run all examples (except `unrestricted.jl`, which will generate a lot of output):
```
$ cd examples
$ ./run_all.sh
```

## Output

- The following output files will be generated when you run each example:
* ${example}_${solver}_optmap.png: the dense map
* ${example}_${solver}_optmap_X.csv: the high-dimensional coordinates of the points used to create the map
* ${example}_${solver}_optmap_C.csv: constraint feasibility for each point
* ${example}_${solver}_optmap_Z.csv: value of the high-dimensional function evaluated at each point in X, in square matrix format (one element per pixel)
* ${example}_${solver}.log: log file from the solver

The output files from `unrestricted.jl` are named slightly differently, but the basic structure and meaning of the contents are the same.
