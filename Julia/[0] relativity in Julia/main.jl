using Plots
using DifferentialEquations
using DiffEqDiffTools

# import/include our files using paths relative to the directory "main.jl"
include(joinpath(@__DIR__, "scenarios", "lambda_L_L_M", "variables.jl"))
include(joinpath(@__DIR__, "scenarios", "lambda_L_L_M", "f_functions.jl"))
include(joinpath(@__DIR__, "scenarios", "lambda_L_L_M", "scenario_data.jl"))
include(joinpath(@__DIR__, "shared", "derivative_functions.jl"))
include(joinpath(@__DIR__, "shared", "potential_extrema.jl"))

# did we correctly import the modules?
using .Variables
using .FFunctions
using .ScenarioData
using .DerivativeFunctionsModule
using .PotentialExtrema

# this is how we extract the "scenario data"
r0   = Float64.(ScenarioData.scene_dict["r0"])  # just to confirm that all initial conditions are of type "Float64"
R0   = ScenarioData.scene_dict["R0"]
mHS  = ScenarioData.scene_dict["mHS"]
rho  = ScenarioData.scene_dict["rho"]
xrho = ScenarioData.scene_dict["xrho"]
T    = ScenarioData.scene_dict["T"]

f    = ScenarioData.scene_dict["f"]
f1   = ScenarioData.scene_dict["f1"]
f2   = ScenarioData.scene_dict["f2"]
f3   = ScenarioData.scene_dict["f3"]
f32  = ScenarioData.scene_dict["f32"]

# this is our instance of the "DerivativeFunctions" structure.
df = DerivativeFunctionsModule.DerivativeFunctions(f, f1, f2, f3, f32, rho, xrho, T)
# this is our instance of the "PotentialDerivative" structure.
pd = PotentialExtrema.PotentialDerivative(f, f1)

# and to plot the potential's derivative, we define
function dVdR_func(R)
    return PotentialExtrema.dVdR(pd, R)
end

# now we can plot the previous function as
# plot(dVdR_func, -1, 10, label="Potential dV/dR") |> display  # we plot within this interval because we already found that for all cases, the "main behaviour" lies here

# then we look for roots of the potential's derivative.
# this will tell us what value "R" (the solution of "R2") approaches
# as the independent parameter "x" (the position) tends to infinity
found_roots = PotentialExtrema.find_roots(dVdR_func, -5, 15, 500)
println("Found extrema at R =", found_roots)

# let's just write a function with no parameters (to be executed automatically when called)
# for the solution of the ODE

#=
function solve_ivp()


    # the function "F!" is the recommended parameter for "ODEProblem"
    prob = ODEProblem(DerivativeFunctionsModule.F!, r0, (1e-6, 1e7), df)
    # maybe use "PFRK87" like in our homework. but if not, try using "Rodas5" (its documentation says it's used for stiff problems)
    sol = solve(prob, Rodas5(), adaptive=false, dt = 0.1, abstol=1e-8, reltol=1e-6)  # as with problems from our homework, the first four parameters are the same. here, we add the last two parameters for better accuracy (as with the working Python code); once we fix the problem :(
    
    #this modification looks more like what we've been using in our homework
    sol_dep = [u[3] for u in sol.u]
    plot(sol.t, sol_dep, xlabel="distance", ylabel="R", title="Integration")

    #plot(sol.t, sol[3, :], xlabel="distance", ylabel="R", title="Integration")
    #display(plot!)
    #println("Integration reached R(", sol.t[end], ") =", sol[3, end])
end

solve_ivp()

=#

function solve_ivp()
    # 1) Wrap F! in an ODEFunction with FD Jacobian
    ode_f = ODEFunction(DerivativeFunctionsModule.F!, jac=DiffEqDiffTools.finite_difference_jacobian)

    # 3) Build and solve with a stiff Rosenbrock method, disabling AD
    prob = ODEProblem(ode_f, r0, (1e-6, 1e7), df)
    sol  = solve(prob,
                 Rodas5(autodiff=false),    # force finite-difference Jacobian
                 abstol=1e-8,
                 reltol=1e-6)

    # 4) Plot R(r) = sol.u[i][3] versus sol.t
    sol_R = getindex.(sol.u, 3)
    plot(sol.t, sol_R;
         xlabel="r",
         ylabel="R(r)",
         title="f(R) integration",
         legend=false)
end

solve_ivp()

# Print test variable from Variables module to confirm proper import.
println(Variables.test_variable)
