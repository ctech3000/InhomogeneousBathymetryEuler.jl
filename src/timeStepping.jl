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

function compute_RHS(g::Vector{Float64},M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, phi_curr::Vector{Float64}, phi_old::Vector{Float64}, phi_oldold::Vector{Float64}, Dt::Real, timeMethod::BackwardDiff)
    return GRAV*g + (2/Dt*M_T1 + 5/Dt^2*M_T0)*phi_curr - (1/(2*Dt)*M_T1 + 4/Dt^2*M_T0)*phi_old + 1/Dt^2*M_T0*phi_oldold
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

function solve_all_timesteps(LHS_matrix::SparseMatrixCSC, LHS_matrix_init::SparseMatrixCSC, M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, domain::AbstractDomain, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64}, time_vec::Vector{Float64}, timeMethod::AbstractTimeSteppingMethod, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, outflow::OutflowBC, D_inflow_boundary::Vector{Vector{Float64}}; save_phi::Union{Bool,Tuple{Int64, Int64, Int64}}=false)
    nχ = length(χs) - 1
    nσ = length(σs) - 1
    Dσ = σs[2] - σs[1]
    nt = length(time_vec)
    Dt = time_vec[2] - time_vec[1]
    phi_old = zeros(Float64, ndofs(dh))       # initial values
    phi_oldold = zeros(Float64, ndofs(dh))    
    phi_curr = zeros(Float64, ndofs(dh))  
    all_etas = Vector{Vector{Float64}}(undef,nt)
    all_etas[1] = zeros(Float64,nχ+1)
    if save_phi !== false
        if isa(save_phi,Bool)
            all_phis = Vector{Matrix{Float64}}(undef,nt)
            all_phis[1] = zeros(Float64,((nχ+1),(nσ+1)))
        else
            skip_χ, skip_σ, skip_t = save_phi 
            all_etas = Vector{Vector{Float64}}(undef,round(Integer,(nt-1)/skip_t+1))
            all_etas[1] = zeros(Float64,round(Integer,nχ/skip_χ+1))
            all_phis = Vector{Matrix{Float64}}(undef,round(Integer,(nt-1)/skip_t+1))
            all_phis[1] = zeros(Float64,((round(Integer,nχ/skip_χ+1)),(round(Integer,nσ/skip_σ+1))))
        end
    end
        
    LHS_lu = lu(LHS_matrix)
    for (t_idx,t_p) in ProgressBar(enumerate(time_vec[2:end]))
        g = assemble_g_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,t_p)
        RHS = compute_RHS(g,M_T0,M_T1,phi_curr,phi_old,phi_oldold,Dt,timeMethod)
        if outflow.type == "Dirichlet"
            apply_dirichlet!(RHS,LHS_matrix_init,zeros(Float64,nσ+1),ch)
        end
        phi_new = LHS_lu\RHS
        eta_new = compute_eta(phi_new,phi_curr,phi_old,domain,trans,Dt,nχ,χs,dh)

        if isa(save_phi,Bool)
            all_etas[t_idx+1] = eta_new
            if save_phi
                phi_nodes = evaluate_at_grid_nodes(dh,phi_new,:phi)
                all_phis[t_idx+1] = reshape(phi_nodes,(nχ+1,nσ+1))
            end
        else
            if round(t_idx/skip_t) == t_idx/skip_t
                phi_nodes_mat = reshape(evaluate_at_grid_nodes(dh,phi_new,:phi),(nχ+1,nσ+1))
                all_phis[round(Integer,t_idx/skip_t)+1] = phi_nodes_mat[1:skip_χ:end,1:skip_σ:end]
                all_etas[round(Integer,t_idx/skip_t)+1] = eta_new[1:skip_χ:end]
            end
        end

        phi_oldold = phi_old
        phi_old = phi_curr
        phi_curr = phi_new
    end

    if isa(save_phi,Tuple{Int64, Int64, Int64}) || save_phi
        return all_etas, all_phis
    else
        return all_etas
    end
end;

function solve_all_timesteps!(sensors::Sensors,LHS_matrix::SparseMatrixCSC, LHS_matrix_init::SparseMatrixCSC, M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, domain::AbstractDomain, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64}, time_vec::Vector{Float64}, timeMethod::AbstractTimeSteppingMethod, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, outflow::OutflowBC, D_inflow_boundary::Vector{Vector{Float64}}; save_phi::Union{Bool,Tuple{Int64, Int64, Int64}}=false)
    if isa(save_phi,Tuple{Int64, Int64, Int64}) || save_phi
        all_etas, all_phis = solve_all_timesteps(LHS_matrix, LHS_matrix_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=save_phi);
        extractSensorData!(sensors,all_etas,χs,save_phi=save_phi)
        return all_etas, all_phis, sensors
    else
        all_etas = solve_all_timesteps(LHS_matrix, LHS_matrix_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=save_phi);
        extractSensorData!(sensors,all_etas,χs,save_phi=save_phi)
        return all_etas, sensors
    end 
end;
