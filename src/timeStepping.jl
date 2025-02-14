#= using Ferrite
using ProgressBars
include("inputParameters.jl")
include("transformations.jl")
include("analyticPotential.jl")
include("assemble.jl")
 =#

function compute_LHS_matrix(K::SparseMatrixCSC, M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, M_T2::SparseMatrixCSC, Dt::Real, timeMethod::BackwardDiff)
    return GRAV*K + M_T2 + 3/(2*Dt)*M_T1 + 2/Dt^2*M_T0
end

function compute_RHS(f::Vector{Float64},M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, phi_curr::Vector{Float64}, phi_old::Vector{Float64}, phi_oldold::Vector{Float64}, Dt::Real, timeMethod::BackwardDiff)
    return GRAV*f + (2/Dt*M_T1 + 5/Dt^2*M_T0)*phi_curr - (1/(2*Dt)*M_T1 + 4/Dt^2*M_T0)*phi_old + 1/Dt^2*M_T0*phi_oldold
end

function compute_σ_derivative_on_free_surface(phi::Vector{Float64},Dσ::Real, nχ::Int)
    phi_dσ = 1/Dσ*(-3/2*phi[1:nχ+1] .+ 2*phi[(2-1)*(nχ+1)+1:2*(nχ+1)] .+ -1/2*phi[(3-1)*(nχ+1)+1:3*(nχ+1)])
    return phi_dσ
end

function compute_new_dirichlet_data(phi_curr::Vector{Float64}, phi_old::Vector{Float64}, domain::DampedDomainProperties, trans::σTransform, χs::Vector{Float64}, Dt::Real,Dσ::Real,nχ::Int)
    phi_surface_old = phi_old[1:nχ+1]
    phi_surface_curr = phi_curr[1:nχ+1]
    phi_dσ = compute_σ_derivative_on_free_surface(phi_curr,Dσ,nχ)
    μ_D_disc = domain.μ_D.(trans.x.(χs))
    bath_disc = eval_bath(trans.tBath,χs)
    phi_surface_new = 1 ./(1/Dt^2 .+ 2*μ_D_disc/Dt).*(phi_surface_curr.*(2/Dt^2 .+ 2*μ_D_disc/Dt .- μ_D_disc.^2) .- phi_surface_old/Dt^2 .- GRAV./bath_disc.*phi_dσ)

    return phi_surface_new
end

function compute_new_dirichlet_data(phi_curr::Vector{Float64}, phi_old::Vector{Float64}, domain::DomainProperties, trans::σTransform, χs::Vector{Float64}, Dt::Real,Dσ::Real,nχ::Int)
    phi_surface_old = phi_old[1:nχ+1]
    phi_surface_curr = phi_curr[1:nχ+1]
    phi_dσ = compute_σ_derivative_on_free_surface(phi_curr,Dσ,nχ)
    bath_disc = eval_bath(trans.tBath,χs)

    phi_surface_new = 2*phi_surface_curr .- phi_surface_old .- Dt^2*GRAV./bath_disc.*phi_dσ

    return phi_surface_new
end

function compute_eta(phi_new::Vector{Float64}, phi_curr::Vector{Float64}, phi_old::Vector{Float64}, domain::DampedDomainProperties, trans::σTransform,Dt::Real,nχ::Int,χs::Vector{Float64}, dh::DofHandler)
    phi_nodes_old = evaluate_at_grid_nodes(dh,phi_old,:phi)
    phi_nodes_curr = evaluate_at_grid_nodes(dh,phi_curr,:phi)
    phi_nodes_new = evaluate_at_grid_nodes(dh,phi_new,:phi)
    phi_surface_old = phi_nodes_old[1:nχ+1]
    phi_surface_curr = phi_nodes_curr[1:nχ+1]
    phi_surface_new = phi_nodes_new[1:nχ+1]
    μ_D_disc = domain.μ_D.(trans.x.(χs))

    eta_new = -1/GRAV*((3/2*phi_surface_new .- 2*phi_surface_curr .+ 1/2*phi_surface_old)/Dt .+ μ_D_disc.*phi_surface_new)

    return eta_new
end

function compute_eta(phi_new::Vector{Float64}, phi_curr::Vector{Float64}, phi_old::Vector{Float64}, domain::DomainProperties, trans::σTransform,Dt::Real,nχ::Int,χs::Vector{Float64}, dh::DofHandler)
    phi_nodes_old = evaluate_at_grid_nodes(dh,phi_old,:phi)
    phi_nodes_curr = evaluate_at_grid_nodes(dh,phi_curr,:phi)
    phi_nodes_new = evaluate_at_grid_nodes(dh,phi_new,:phi)
    phi_surface_old = phi_nodes_old[1:nχ+1]
    phi_surface_curr = phi_nodes_curr[1:nχ+1]
    phi_surface_new = phi_nodes_new[1:nχ+1]

    eta_new = -1/GRAV*((3/2*phi_surface_new .- 2*phi_surface_curr .+ 1/2*phi_surface_old)/Dt)

    return eta_new
end

function solve_all_timesteps(LHS_matrix::SparseMatrixCSC, LHS_matrix_init::SparseMatrixCSC, M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, domain::AbstractDomain, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64}, time_vec::Vector{Float64}, timeMethod::AbstractTimeSteppingMethod, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, outflow::OutflowBC, D_inflow_boundary::Vector{Vector{Float64}}; save_phi::Bool=false)
    nχ = length(χs) - 1
    nσ = length(σs) - 1
    Dσ = σs[2] - σs[1]
    nt = length(time_vec)
    Dt = time_vec[2] - time_vec[1]
    phi_old = zeros(Float64, ((nχ+1)*(nσ+1)))       # initial values
    phi_oldold = zeros(Float64, ((nχ+1)*(nσ+1)))    
    phi_curr = zeros(Float64, ((nχ+1)*(nσ+1)))  
    all_etas = Vector{Vector{Float64}}(undef,nt)
    all_etas[1] = zeros(Float64,nχ+1)
    if save_phi
        all_phis = Vector{Matrix{Float64}}(undef,nt)
        all_phis[1] = zeros(Float64,((nχ+1),(nσ+1)))
    end
        
           
    for (t_idx,t_p) in ProgressBar(enumerate(time_vec[2:end]))
        f = assemble_f_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,t_p)
        RHS = compute_RHS(f,M_T0,M_T1,phi_curr,phi_old,phi_oldold,Dt,timeMethod)
        if outflow.type == "Dirichlet"
            apply_dirichlet!(RHS,LHS_matrix_init,zeros(Float64,nσ+1),ch)
        end
        phi_new = LHS_matrix\RHS
        eta_new = compute_eta(phi_new,phi_curr,phi_old,domain,trans,Dt,nχ,χs,dh)

        all_etas[t_idx+1] = eta_new
        if save_phi
            phi_nodes = evaluate_at_grid_nodes(dh,phi_new,:phi)
            all_phis[t_idx+1] = reshape(phi_nodes,(nχ+1,nσ+1))
        end

        phi_oldold = phi_old
        phi_old = phi_curr
        phi_curr = phi_new
    end

    if save_phi
        return all_etas, all_phis
    else
        return all_etas
    end
end;

#= function solve_all_timesteps(K::SparseMatrixCSC, K_init::SparseMatrixCSC, domain::AbstractDomain, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64}, time_vec::Vector{Float64}, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, D_inflow_boundary::Vector{Vector{Float64}}; save_phi::Bool=false)
    nχ = length(χs) - 1
    nσ = length(σs) - 1
    Dσ = σs[2] - σs[1]
    nt = length(time_vec)
    Dt = time_vec[2] - time_vec[1]
    all_etas = Vector{Vector{Float64}}(undef,nt)
    if save_phi
        all_phis = Vector{Matrix{Float64}}(undef,nt)
        all_phi_surface = Vector{Vector{Float64}}(undef,nt)
    end
    phi_old = zeros(Float64, ((nχ+1)*(nσ+1)))       # this is an initial value
    phi_curr = zeros(Float64, ((nχ+1)*(nσ+1)))      # this is not an initial value, it will be overwritten immediately 
    phi_surface_curr = zeros(Float64,nχ+1)          # this is an initial value   
    


    for (t_idx,t_p) in ProgressBar(enumerate(time_vec))
        f = assemble_f_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,t_p)
        apply_dirichlet!(f,K_init,phi_surface_curr,ch)
        phi_curr = evaluate_at_grid_nodes(dh,K\f,:phi)
        eta_curr = compute_eta(phi_old,phi_curr,domain,trans,Dt,nχ,χs::Vector{Float64})

        phi_surface_new = compute_new_dirichlet_data(phi_curr,phi_old,domain,trans,χs,Dt,Dσ,nχ)
        
        all_etas[t_idx] = eta_curr
        if save_phi
            all_phis[t_idx] = reshape(phi_curr,(nχ+1,nσ+1))
            all_phi_surface[t_idx] = phi_surface_new
        end
        phi_surface_curr = phi_surface_new
        phi_old = phi_curr
    end

    if save_phi
        return all_etas, all_phis, all_phi_surface
    else
        return all_etas
    end
end

function solve_one_timestep(K::SparseMatrixCSC, K_init::SparseMatrixCSC, dirichlet_data::Vector{Float64}, domain::AbstractDomain, trans::σTransform, t_p::Float64, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, D_inflow_boundary::Vector{Vector{Float64}})
    f = assemble_f_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,t_p)
    apply_dirichlet!(f,K_init,dirichlet_data,ch)
    phi_curr = evaluate_at_grid_nodes(dh,K\f,:phi)
    return phi_curr
end =#