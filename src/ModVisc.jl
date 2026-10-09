module ModVisc

export AbstractGeometry,
       Film,
       Tube,
       Spheroid,
       MeltFraction,
       PoreFraction,
       ViscousModel,
       ElasticModel,
       SolverOptions,
       ViscosityResult,
       ElasticResult,
       viscosity,
       viscosity!,
       viscosity_profile,
       viscosity_profile!,
       elasticity!,
       elasticity_profile,
       elasticity_profile!,
       poisson_ratio,
       elasticity,
       spheroid_fit,
       spheroid_fit!

include("geometries.jl")
include("parameters.jl")
include("kernels.jl")
include("viscous.jl")
include("elastic.jl")
include("curves.jl")

end
