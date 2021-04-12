function get_cigar_funcs(dimensions)
    fs = Dict()

    function cigar()
        f(x) = x[1]^2 + 10^6 * sum(x[2:end].^2)
        return f
    end

    for d in dimensions
        fs[d] = cigar()
    end

    return fs
end

function get_linear_funcs(dimensions)
    fs = Dict()

    function linear()
        f(x) = sum(x)
        return f
    end

    for d in dimensions
        fs[d] = linear()
    end

    return fs
end

function get_box_funcs(dimensions, alpha=1/2)
    fs = Dict()

    function box()
        f(x) = sum(abs.(x)) <= alpha
        return f
    end

    for d in dimensions
        fs[d] = box()
    end

    return fs
end

function get_step_funcs(dimensions)
    fs = Dict()

    function step()
        f(x) = begin
            v = sum(x)
            if v < 0 
                return 0
            elseif v < 2
                return 2
            elseif v < 4
                return 4
            else
                return 5
            end
        end
           
        return f
    end

    for d in dimensions
        fs[d] = step()
    end

    return fs
end

function get_sinc_funcs(dimensions)
    fs = Dict()

    function msinc()
        f(x) = *(sinc.(x)...)
        return f
    end

    for d in dimensions
        fs[d] = msinc()
    end

    return fs
end

function get_sphere_funcs(dimensions)
    fs = Dict()

    function sphere()
        f(x) = sum(x.^2)
        return f
    end

    for d in dimensions
        fs[d] = sphere()
    end

    return fs
end

function get_rastrigin_funcs(dimensions, A=10.0)
    fs = Dict()

    function rastrigrin(n)
        f(x) = A*n + sum(x.^2 - A * cos.(2π * x))
        return f
    end

    for d in dimensions
        fs[d] = rastrigrin(d)
    end

    return fs
end

function get_rosenbrock_funcs(dimensions)
    fs = Dict()

    function rosenbrock(n)
        f(x) = sum([100 * (x[i+1] - x[i]^2) + (1 - x[i])^2 for i in 1:(n-1)])
        return f
    end

    for d in dimensions
        fs[d] = rosenbrock(d)
    end

    return fs
end


function get_styblinskitang_funcs(dimensions)
    fs = Dict()

    function styblinskitang()
        f(x) = sum(x.^4 - 16x.^2 + 5x) / 2
        return f
    end

    for d in dimensions
        fs[d] = styblinskitang()
    end

    return fs
end
