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
    if x == 0
        return 0.0
    elseif P > 0
        return (n / (x * (2 * df.f1(R) + x * R1 * df.f2(R)))) *
               ((x^2) * m * (df.f(R) - R * df.f1(R) + 2 * 8 * π * P) +
                2 * (m - 1) * df.f1(R) - 4 * x * R1 * df.f2(R))
    else
        return (n / (x * (2 * df.f1(R) + x * R1 * df.f2(R)))) *
               ((x^2) * m * (df.f(R) - R * df.f1(R)) +
                2 * (m - 1) * df.f1(R) - 4 * x * R1 * df.f2(R))
    end
end

function m1(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    if x == 0
        return 0.0
    elseif P > 0
        return (m / (x * (2 * df.f1(R) + x * R1 * df.f2(R)))) *
               (2 * df.f1(R) * (1 - m) + 16 * π * m * (x^2) * df.rho +
                ((m * (x^2)) / 3) * (R * df.f1(R) + df.f(R) + 16 * π * df.T) +
                x * R1 * (df.f2(R) / df.f1(R)) *
                (((m * (x^2)) / 3) * (2 * R * df.f1(R) - df.f(R) + 8 * π * df.T) -
                 8 * π * m * (x^2) * (-df.rho + P) +
                 2 * (1 - m) * df.f1(R) + 2 * x * R1 * df.f2(R)))
    else
        return (m / (x * (2 * df.f1(R) + x * R1 * df.f2(R)))) *
               (2 * df.f1(R) * (1 - m) +
                ((m * (x^2)) / 3) * (R * df.f1(R) + df.f(R)) +
                x * R1 * (df.f2(R) / df.f1(R)) *
                (((m * (x^2)) / 3) * (2 * R * df.f1(R) - df.f(R)) +
                 2 * (1 - m) * df.f1(R) + 2 * x * R1 * df.f2(R)))
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
    if x == 0
        return (2 * df.f(R) - df.f1(R) * R + 8 * π * df.T) / (9 * df.f2(R))
    elseif P > 0
        m1_val = m1(df, x, n, m, R, R1, P, Ms, Mb)
        n1_val = n1(df, x, n, m, R, R1, P, Ms, Mb)
        return (m * (8 * π * df.T + 2 * df.f(R) - df.f1(R) * R)) / (3 * df.f2(R)) -
               df.f32(R) * (R1^2) + m1_val * (R1 / (2 * m)) -
               ((n1_val / (2 * n)) + (2 / x)) * R1
    elseif x <= 0.1 && P <= 0
        m1_val = m1(df, x, n, m, R, R1, P, Ms, Mb)
        n1_val = n1(df, x, n, m, R, R1, P, Ms, Mb)
        return (m * (2 * df.f(R) - df.f1(R) * R)) / (3 * df.f2(R)) -
               df.f32(R) * (R1^2) + m1_val * (R1 / (2 * m)) -
               ((n1_val / (2 * n)) + (2 / x)) * R1
    else
        return 0.0
    end
end

function DP(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    if x == 0
        return 0.0
    elseif P > 0
        n1_val = n1(df, x, n, m, R, R1, P, Ms, Mb)
        return - (df.rho + P) * (n1_val / (2 * n))
    else
        return 0.0
    end
end

function DMs(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    if x == 0
        return 0.0
    elseif P > 0
        return 4 * π * df.rho * (x^2)
    else
        return 0.0
    end
end

function DMb(df::DerivativeFunctions, x, n, m, R, R1, P, Ms, Mb)
    if x == 0
        return 0.0
    elseif P > 0
        return 4 * π * df.rho * (x^2) * df.xrho
    else
        return 0.0
    end
end

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

function G(df::DerivativeFunctions, u, r)
    x = u / (1 - u)
    factor = (1 - u)^(-2)
    return factor * F(df, x, r)
end

end
