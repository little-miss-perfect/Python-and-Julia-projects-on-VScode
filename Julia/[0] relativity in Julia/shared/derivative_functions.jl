module DerivativeFunctionsModule

export DerivativeFunctions, n1, m1, DR, R2, DP, DMs, DMb, F, F!, G

struct DerivativeFunctions
    f::Function
    f1::Function
    f2::Function
    f3::Function
    f32::Function
    rho::Float64
    xrho::Float64
    T::Float64
end

function n1(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    
    coeff = n / (x * (2 * df.f1(R) + x * R1 * df.f2(R)))
    sum1 = m * (x^2) * (df.f(R) - R * df.f1(R))
    sum2 = 2 * df.f1(R) * (m - 1) - 4 * x * R1 * df.f2(R)
    kappa = 8 * π
    
    if x == 0
        return 0.0
    elseif P > 0  # this case is what's presented in the article
        return coeff *
               (sum1 + 2 * kappa * P) +
                sum2
    else  # here "P=0"
        return coeff *
               sum1 +
                sum2
    end
end

function m1(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)

    coeff1 = m / (x * (2 * df.f1(R) + x * R1 * df.f2(R)))
    coeff2 = x * R1 * (df.f2(R) / df.f1(R))
    coeff3 = (m * (x^2)) / 3
    kappa = 8 * π

    sum1 = 2 * df.f1(R) * (1 - m)
    sum2 = R * df.f1(R) + df.f(R)
    sum3 = 2 * R * df.f1(R) - df.f(R)
    sum4 = 2 * x * R1 * df.f2(R)

    if x == 0
        return 0.0
    elseif P > 0  # this case is what's presented in the article
        return coeff1 *
               (sum1 - 2 * kappa * m * (x^2) * (-df.rho) +
                coeff3 * (sum2 + 2 * kappa * df.T) +
                coeff2 *
                (coeff3 * (sum3 + kappa * df.T) -
                 kappa * m * (x^2) * (-df.rho + P) +
                 sum1 + sum4))
    else  # here "P=0"
        return coeff1 *
               (sum1 +
                coeff3 * sum2 +
                coeff2 *
                (coeff3 * sum3 +
                sum1 + sum4))
    end
end

function DR(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    if x == 0
        return 0.0
    elseif P > 0 || (x <= 0.1 && P <= 0)
        return R1
    else
        return 0.0
    end
end

function R2(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)

    m1_val = m1(df, x, n, m, R, R1, P, Ms, Mb)
    n1_val = n1(df, x, n, m, R, R1, P, Ms, Mb)

    kappa = 8 * π

    sum1 = 2 * df.f(R) - R * df.f1(R)
    sum2 = 3 * df.f3(R) * (R1^2)
    sum3 = (m1_val/(2*m)) - (n1_val/(2*n)) - 2/x

    coeff1 = 1 / (3 * df.f2(R))
    
    if x == 0  # defined on page "3" of the article
        return (sum1 + kappa * df.T) / (9 * df.f2(R))

    elseif P > 0  # this case is what's presented in the article
        return coeff1 * (m * (kappa * df.T + sum1) - sum2) +
               sum3 * R1

    elseif x <= 0.1 && P <= 0  # here "T=0"
        m1_val = m1(df, x, n, m, R, R1, P, Ms, Mb)
        n1_val = n1(df, x, n, m, R, R1, P, Ms, Mb)
        return coeff1 * (m * (sum1)) -
               df.f32(R) * (R1^2) + 
               sum3 * R1

    else
        return 0.0
    end
end

function DP(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    if x == 0  # by how we define "DP", on page "3" of the article, if "x=0", then "n1=0"
        return 0.0  # since "n1=0", and "DP" is defined as a product with "n1", then "DP=0"
    elseif P > 0  # defined on page "3" of the article
        n1_val = n1(df, x, n, m, R, R1, P, Ms, Mb)
        return - (df.rho + P) * (n1_val / (2 * n))
    else
        return 0.0
    end
end

function DMs(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)

    kappa = 8 * π

    if x == 0
        return 0.0
    elseif P > 0
        return (1/2) * kappa * (x^2) * df.rho  # this definition is given in the original Fortran program
    else
        return 0.0
    end
end

function DMb(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)    
    return DMs(df, x, n, m, R, R1, P, Ms, Mb) * df.xrho  # this definition is given in the original Fortran program
end

function F!(du, u, p, x)
    # unpack the state vector once
    n, m, R, R1, P, Ms, Mb = u

    # write each derivative in-place
    du[1] = n1(p, x, n, m, R, R1, P, Ms, Mb)
    du[2] = m1(p, x, n, m, R, R1, P, Ms, Mb)
    du[3] = DR(p, x, n, m, R, R1, P, Ms, Mb)
    du[4] = R2(p, x, n, m, R, R1, P, Ms, Mb)
    du[5] = DP(p, x, n, m, R, R1, P, Ms, Mb)
    du[6] = DMs(p, x, n, m, R, R1, P, Ms, Mb)
    du[7] = DMb(p, x, n, m, R, R1, P, Ms, Mb)

    return nothing
end

#=
function F(df::DerivativeFunctions, x, r)
    n   = r[1]
    m   = r[2]
    R   = r[3]
    R1  = r[4]
    P   = r[5]
    Ms  = r[6]
    Mb  = r[7]
    Fvec = Vector{Float64}(undef, 7)
    Fvec[1] = n1(df, x, n, m, R, R1, P, Ms, Mb)
    Fvec[2] = m1(df, x, n, m, R, R1, P, Ms, Mb)
    Fvec[3] = DR(df, x, n, m, R, R1, P, Ms, Mb)
    Fvec[4] = R2(df, x, n, m, R, R1, P, Ms, Mb)
    Fvec[5] = DP(df, x, n, m, R, R1, P, Ms, Mb)
    Fvec[6] = DMs(df, x, n, m, R, R1, P, Ms, Mb)
    Fvec[7] = DMb(df, x, n, m, R, R1, P, Ms, Mb)
    return Fvec
end

function F!(du, u, p, t)
    deriv = F(p, t, u)
    for i in eachindex(du)
        du[i] = deriv[i]
    end
    return nothing
end

function G!(df::DerivativeFunctions, u, r)
    x = u / (1 - u)
    factor = (1 - u)^(-2)
    return factor * F(df, x, r)
end
=#

end  # to end the module
