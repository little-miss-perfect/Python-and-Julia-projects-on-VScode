module PotentialExtrema

using Roots

export PotentialDerivative, dVdR, find_roots

struct PotentialDerivative
    f::Function
    f1::Function
end

function dVdR(pd::PotentialDerivative, R)
    """
    dV/dR = (1/3) * (2*f(R) - R*f1(R))
    """
    return (1/3) * (2 * pd.f(R) - pd.f1(R) * R)
end

function find_roots(f, x_start, x_end, num_points=500)
    x_grid = range(x_start, x_end, length=num_points)
    roots = Float64[]
    for i in 1:length(x_grid)-1
        x_left = x_grid[i]
        x_right = x_grid[i+1]
        f_left = f(x_left)
        f_right = f(x_right)
        tol = 1e-12
        if abs(f_left) < tol
            push!(roots, x_left)
            continue
        end
        if abs(f_right) < tol
            push!(roots, x_right)
            continue
        end
        if f_left * f_right < 0
            try
                root = find_zero(f, (x_left, x_right), Bisection())
                push!(roots, root)
            catch e
                # If no root is found, continue silently.
            end
        end
    end
    return roots
end

end  # module PotentialExtrema
