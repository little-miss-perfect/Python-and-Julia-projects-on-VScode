module ScenarioData

using ..Variables
using ..FFunctions

export scene_dict

const scene_dict = Dict(
    "r0"   => Variables.r0,
    "R0"   => Variables.R0,
    "mHS"  => Variables.mHS,
    "rho"  => Variables.rho,
    "xrho" => Variables.xrho,
    "T"    => Variables.T,
    "f"    => FFunctions.f,
    "f1"   => FFunctions.f1,
    "f2"   => FFunctions.f2,
    "f3"   => FFunctions.f3,
    "f32"  => FFunctions.f32
)

end  # module ScenarioData
