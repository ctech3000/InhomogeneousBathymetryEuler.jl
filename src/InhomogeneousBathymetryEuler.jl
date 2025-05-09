module InhomogeneousBathymetryEuler

using Interpolations, Ferrite, SparseArrays, ProgressBars, JLD2, BSplineKit, PlutoUI, FFTW, DelimitedFiles, LinearAlgebra, GLMakie, HDF5, LaTeXStrings

export computeWavenumber, analyticPotential, analyticPotential_dx, analyticPotential_dz, transformedAnalyticPotential, transformedAnalyticPotential_dχ, transformedAnalyticPotential_dσ, analyticPotential_dt, transformedAnalyticPotential_dt, analyticEtaSimulated
export dirichlet_from_discretized_data, assemble_K_element!, assemble_K_global, assemble_Ms_global, assemble_g_element!, insert_into_f!, assemble_g_global, apply_dirichlet!, meandiag, init_K_M
export get_boundary_coordinates, discretizeTransformedDomain
export firstDerivative, dofToCoordinate
export AbstractWave, SimpleWave, IrregWave, getFreq, getAmp, getWavenumber, Bathymetry, eval_bath, AbstractDomain, DomainProperties, DampedDomainProperties, RelaxedDampedDomainProperties, BackwardDiff, OutflowBC
export compute_σ_derivative_on_free_surface, compute_new_dirichlet_data, compute_eta, solve_all_timesteps, solve_all_timesteps!, solve_one_timestep, computePhysDomainMask, computeEnergy
export σTransform, compute_B_D
export computeDerivativeOnBoundary, laplace
export GRAV
export TrueSolution, assemble_f_global, assemble_g_global, coefficientVector, assemble_h_global, computeBathNormal, assemble_l_global
export Sensors, extractSensorData!
export computeError, computeErrorL2, computeErrorMax, computeEOC
export plotSurfaceOverTime, plotSurfaceOverTime!, plotSensorData, plotSensorData!, plotSensorData_alt, loadSWE
#export plotSensorData, plotSensorData!

include("helpfulFunctions.jl")
include("inputParameters.jl")
include("transformations.jl")
include("sensors.jl")
include("analyticPotential.jl")
include("domainDiscretization.jl")
include("assemble.jl")
include("timeStepping.jl")
include("testSolution.jl")
include("manufacturedSolution.jl")
include("errorComputation.jl")
include("plotFunctions.jl")

end
