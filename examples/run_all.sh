#!/bin/bash

for test_script in diet.jl prod_sched.jl knapsack.jl nonlinear_sphere.jl nonlinear_styblinski_tang.jl
do
    echo "######################################################"
    echo "$test_script"
    echo "######################################################"
    julia -J ../sys_optmap.so $test_script
done
