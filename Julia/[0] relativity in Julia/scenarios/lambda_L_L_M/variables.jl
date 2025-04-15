module Variables

using CSV
using DataFrames

# Export the symbols needed in other modules.
export lambda_staro, R_star, beta, this_dir, data_cond_path, data_gas_path,
       data_cond, data_gas, x0, n0, n10, m0, m10, R0, DR0, Ms0, delta,
       rho, nHS, mHS, f1hoy, c2, c1, P0, Mb0, G, c, M_sol, rho_nuc, rho_c, xrho,
       Po, T, r0, f0, f10, f20, f30, test_variable

# --------------------------
# Model Parameters and Paths
# --------------------------
const lambda_staro = 1.56
const R_star = 1.0
const beta = 1

# the absolute path to this script's directory
const this_dir = @__DIR__

# Build paths to the data files (portable)
const data_cond_path = joinpath(this_dir, "data", "init_cond.dat")
const data_gas_path  = joinpath(this_dir, "data", "init_cond.dat")  # Adjust if file name differs

# --------------------------
# Read Data Files into DataFrames
# --------------------------
const data_cond = CSV.File(data_cond_path) |> DataFrame
const data_gas  = CSV.File(data_gas_path) |> DataFrame

# --------------------------
# Assign Variables from "data_cond"
# --------------------------
const x0 = data_cond[1, 1]
const n0 = data_cond[2, 1]
const n10 = data_cond[3, 1]
const m0 = data_cond[4, 1]
const m10 = data_cond[5, 1]
const R0 = data_cond[6, 1]   # This is the value we want to find via a shooting method
const DR0 = data_cond[7, 1]
const Ms0 = data_cond[8, 1]
const delta = data_cond[9, 1]

# --------------------------
# Assign Variables from "data_gas"
# --------------------------
const rho = data_gas[1, 1]
const nHS = data_gas[2, 1]
const mHS = data_gas[3, 1]  # After squaring, used in f2(R)
const f1hoy = data_gas[4, 1]

# --------------------------
# Define Other Variables
# --------------------------
const c2 = (-1 * nHS * 19 * 41^(-1 * nHS - 1)) / f1hoy
const c1 = 19 * c2
const P0 = 0.1 * rho

const Mb0 = 1.66e-27   # initial baryonic mass?
const G = 1
const c = 1
const M_sol = 1

const rho_nuc = 1.66e17
const rho_c = (c^8) / (M_sol^2 * G^3)
const xrho = (rho_nuc * c^2) / rho_c

const Po = 0.1 * rho
const T = -rho + 3 * Po

# --------------------------
# Initial Conditions Vector
# --------------------------
const r0 = [n0, m0, R0, DR0, P0, Ms0, Mb0]

# --------------------------
# Define Functions for Numerical Method Comparison
# --------------------------
const f0 = R0 - (c1 * mHS^2 * (R0 / mHS^2)^nHS) / (1 + c2 * (R0 / mHS^2)^nHS)
const f10 = 1 - (c1 * nHS * (R0 / mHS^2)^(-1 + nHS)) / ((1 + c2 * (R0 / mHS^2)^nHS)^2)
const f20 = (c1 * mHS^2 * nHS * (R0 / mHS^2)^nHS * (1 + c2 * (R0 / mHS^2)^nHS + nHS * (-1 + c2 * (R0 / mHS^2)^nHS))) /
            (R0^2 * (1 + c2 * (R0 / mHS^2)^nHS)^3)
const f30 = -(c1 * mHS^2 * nHS * (R0 / mHS^2)^nHS * (2 * (1 + c2 * (R0 / mHS^2)^nHS)^2 +
            3 * nHS * (-1 + c2^2 * (R0 / mHS^2)^(2 * nHS)) +
            nHS^2 * (1 - 4 * c2 * (R0 / mHS^2)^nHS + c2^2 * (R0 / mHS^2)^(2 * nHS)))) /
            (R0^3 * (1 + c2 * (R0 / mHS^2)^nHS)^4)

# Test variable to check proper import.
const test_variable = "we got the information from the original file! after nesting these files :)"

end  # module Variables
