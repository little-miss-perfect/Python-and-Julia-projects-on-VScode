module FFunctions

using ..Variables  # Use two dots to move up one level

export f, f1, f2, f3, f32

function f(R)
    return R - (Variables.lambda_staro * Variables.R_star) * (1 - (1 + (R / Variables.R_star)^2)^(-Variables.beta))
end

function f1(R)
    """
    First derivative f'(R).
    """
    R_tilde = R / Variables.R_star
    return 1 - 2 * Variables.lambda_staro * Variables.beta * R_tilde * (1 + R_tilde^2)^(-Variables.beta - 1)
end

function f2(R)
    """
    Second derivative f''(R).
    """
    R_tilde = R / Variables.R_star
    factor = -2 * Variables.lambda_staro * Variables.beta / Variables.R_star
    bracket = (1 + R_tilde^2) + 2 * (-Variables.beta - 1) * R * R_tilde / Variables.R_star
    return factor * ((1 + R_tilde^2)^(-Variables.beta - 2) * bracket)
end

function f3(R)
    """
    Third derivative f'''(R).
    """
    R_tilde = R / Variables.R_star
    A = (1 + R_tilde^2)^(-Variables.beta - 2)
    B = (1 + R_tilde^2) + 2 * (-Variables.beta - 1) * R * R_tilde / Variables.R_star
    A_prime = -2 * (Variables.beta + 2) * R_tilde / Variables.R_star * (1 + R_tilde^2)^(-Variables.beta - 3)
    B_prime = -2 * (2 * Variables.beta + 1) * R / (Variables.R_star^2)
    g2 = A_prime * B + A * B_prime
    factor = -(2 * Variables.lambda_staro * Variables.beta / Variables.R_star)
    return factor * g2
end

function f32(R)
    return f3(R) / f2(R)
end

end  # module FFunctions
