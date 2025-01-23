module InhomogeneousBathymetryEuler

using Interpolations, Ferrite, SparseArrays, ProgressBars

export computeWaveNumber, analyticPotential, analyticPotential_dx, analyticPotential_dz, transformedAnalyticPotential, transformedAnalyticPotential_dχ, transformedAnalyticPotential_dσ
export dirichlet_from_discretized_data, assemble_K_element!, assemble_K_global, assemble_f_element!, insert_into_f!, assemble_f_global, init_K_f, apply_dirichlet!, meandiag
export get_boundary_coordinates, discretizeTransformedDomain
export firstDerivative
export AbstractWave, SimpleWave, Bathymetry, eval_bath, AbstractDomain, DomainProperties, DampedDomainProperties
export compute_σ_derivative_on_free_surface, compute_new_dirichlet_data, compute_eta, solve_all_timesteps, solve_one_timestep
export σTransform, compute_B_D
export computeDerivativeOnBoundary, laplace

include("helpfulFunctions.jl")
include("inputParameters.jl")
include("transformations.jl")
include("analyticPotential.jl")
include("domainDiscretization.jl")
include("assemble.jl")
include("timeStepping.jl")
include("testSolution.jl")

end
