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

function relaxation_dofs_coords(dh::Ferrite.DofHandler, trans::σTransform, domain::RelaxedDampedDomainProperties; apply_full_relaxation::Bool=false)
    node_ids = 1:length(dh.grid.nodes)
    RX_dofs = Int[]
    RX_coords = Vector{Float64}[]
    x_RX = domain.x_RX
    x_L = domain.x_L

    for node_id in node_ids
        dof_coord = dh.grid.nodes[node_id].x
        x = trans.x(dof_coord[1])
        z = trans.z(dof_coord...)
        if x <= x_L && (abs(z) < 10^-13 || apply_full_relaxation)
            vertexid = nodeid_to_vertexindex(dh.grid, node_id)
            dof = vertexdofs(dh, vertexid)[1]
            push!(RX_dofs,dof)
            push!(RX_coords,[x,z])
        end
    end
    return RX_dofs, RX_coords
end

function relaxation_dofs_coords(dh::Ferrite.DofHandler, trans::σTransform, domain::AbstractDomain; apply_full_relaxation::Bool=false)
    return [0],[[0.0]]
end

function apply_relaxation(sol_num::Vector{Float64},t::Real,h::Real,RX_dofs::Vector{Int},RX_coords::Vector{Vector{Float64}},wave::AbstractWave,domain::RelaxedDampedDomainProperties)
    relaxed_sol = deepcopy(sol_num)
    for (i,dof) in enumerate(RX_dofs)
        x,z = RX_coords[i]
        sol_ana = analyticPotential(x,z,t,h,wave)
        relaxed_sol[dof] = domain.ramp_RX(x)*sol_num[dof] + (1 - domain.ramp_RX(x))*sol_ana
        #relaxed_sol[dof] = domain.ramp_RX(x)*sol_num[dof] + (1 - domain.ramp_RX(x))*0
    end
    return relaxed_sol
end

function apply_relaxation(sol_num::Vector{Float64},t::Real,h::Real,RX_dofs::Vector{Int},RX_coords::Vector{Vector{Float64}},wave::AbstractWave,domain::AbstractDomain)
    return sol_num
end

function compute_eta(phi_new::Vector{Float64}, phi_curr::Vector{Float64}, phi_old::Vector{Float64}, domain::Union{DampedDomainProperties,RelaxedDampedDomainProperties}, trans::σTransform,Dt::Real,nχ::Int,χs::Vector{Float64}, dh::DofHandler)
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

function find_physical_ind(domain::AbstractDomain,trans::σTransform,χs::Vector{Float64})
    χ_L = trans.χ(domain.x_L)
    if domain isa DomainProperties
        χ_D = trans.χ(domain.x_R)
    else
        χ_D = trans.χ(domain.x_D)
    end
    idx_χL = findall(χ->abs(χ-χ_L)==minimum(abs.(χs.-χ_L)),χs)[1]
    idx_χD = findall(χ->abs(χ-χ_D)==minimum(abs.(χs.-χ_D)),χs)[1]
    return idx_χL, idx_χD
end

function computePhysDomainMask(dh::DofHandler, trans::σTransform, domain::AbstractDomain, Dχ::Real)
    Dx = -trans.x(Dχ)
    if domain isa DomainProperties
        rightBound = domain.x_R
        nextX = domain.x_L-1
    else
        rightBound = domain.x_D 
        nextX = rightBound + Dx
    end
    etaMask = zeros(ndofs(dh))
    phiMask = zeros(ndofs(dh))
    for nodeId = 1:ndofs(dh)
        vertexid = nodeid_to_vertexindex(dh.grid, nodeId)
        dof = vertexdofs(dh, vertexid)[1]
        dof_coord = dh.grid.nodes[nodeId].x
        x = trans.x(dof_coord[1])
        z = trans.z(dof_coord...)
        if x >= domain.x_L && x <= rightBound
            phiMask[dof] = 1
            if abs(z) < 10^-12
                etaMask[dof] = 1
            end
        elseif abs(x-nextX) < 1/2*Dx
            phiMask[dof] = -1
            if abs(z) < 10^-12
                etaMask[dof] = -1
            end
        end
    end
    return phiMask, etaMask
    
end

function computeEnergy(phi::Vector{Float64},eta::Vector{Float64},K::SparseMatrixCSC,M_T0::SparseMatrixCSC)
    return 1/2*transpose(phi)*K*phi + 1/2*GRAV*transpose(eta)*M_T0*eta
end

function solve_all_timesteps(LHS_matrix::SparseMatrixCSC, LHS_matrix_init::SparseMatrixCSC, K_init::SparseMatrixCSC, M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, domain::AbstractDomain, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64}, time_vec::Vector{Float64}, timeMethod::AbstractTimeSteppingMethod, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, outflow::OutflowBC, D_inflow_boundary::Vector{Vector{Float64}}; save_phi::Union{Bool,Tuple{Int64, Int64, Int64}}=false, apply_full_relaxation::Bool=false, energy::Union{Vector{Float64},Nothing}=nothing)
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
    LHS_lu = lu(LHS_matrix)
    RX_dofs, RX_coords = relaxation_dofs_coords(dh,trans,domain,apply_full_relaxation=apply_full_relaxation)

    if save_phi !== false
        if isa(save_phi,Bool)
            all_phis = Vector{Matrix{Float64}}(undef,nt)
            all_phis[1] = zeros(Float64,((nχ+1),(nσ+1)))
        else
            skip_χ, skip_σ, skip_t = save_phi 
            idx_χL, idx_χD = find_physical_ind(domain,trans,χs)
            nχ_phys = idx_χL - idx_χD 
            all_etas = Vector{Vector{Float64}}(undef,round(Integer,(nt-1)/skip_t+1))
            all_etas[1] = zeros(Float64,round(Integer,nχ_phys/skip_χ+1))
            all_phis = Vector{Matrix{Float64}}(undef,round(Integer,(nt-1)/skip_t+1))
            all_phis[1] = zeros(Float64,((round(Integer,nχ_phys/skip_χ+1)),(round(Integer,nσ/skip_σ+1))))
        end
    end
    if !isnothing(energy)
        push!(energy,0.0)
        phiMask, etaMask =  computePhysDomainMask(dh,trans,domain,χs[2]-χs[1])
    end
    
    for (t_idx,t_p) in ProgressBar(enumerate(time_vec[2:end]))
        g = assemble_g_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,t_p)
        RHS = compute_RHS(g,M_T0,M_T1,phi_curr,phi_old,phi_oldold,Dt,timeMethod)
        if outflow.type == "Dirichlet"
            apply_dirichlet!(RHS,LHS_matrix_init,zeros(Float64,nσ+1),ch)
        end
        phi_new = LHS_lu\RHS
        phi_new = apply_relaxation(phi_new,t_p,-domain.b_L,RX_dofs,RX_coords,domain.wave,domain)
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
                all_phis[round(Integer,t_idx/skip_t)+1] = phi_nodes_mat[idx_χD:skip_χ:idx_χL,1:skip_σ:end]
                all_etas[round(Integer,t_idx/skip_t)+1] = eta_new[idx_χD:skip_χ:idx_χL]
            end
        end

        if !isnothing(energy)
            eta_coeff = coefficientVector(dh,eta_new)
            eta_coeff_phys = eta_coeff.*etaMask #restrict phi and eta to (almost only) physical domain
            phi_coeff_phys = phi_new.*phiMask
            energy_new = computeEnergy(phi_coeff_phys,eta_coeff_phys,K_init,M_T0)
            push!(energy,energy_new)
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

function solve_all_timesteps!(sensors::Sensors,LHS_matrix::SparseMatrixCSC, LHS_matrix_init::SparseMatrixCSC, K_init::SparseMatrixCSC, M_T0::SparseMatrixCSC, M_T1::SparseMatrixCSC, domain::AbstractDomain, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64}, time_vec::Vector{Float64}, timeMethod::AbstractTimeSteppingMethod, facetvalues::FacetValues, dh::DofHandler, ch::ConstraintHandler, outflow::OutflowBC, D_inflow_boundary::Vector{Vector{Float64}}; save_phi::Union{Bool,Tuple{Int64, Int64, Int64}}=false, apply_full_relaxation::Bool=false)
    if isa(save_phi,Tuple{Int64, Int64, Int64}) || save_phi
        all_etas, all_phis = solve_all_timesteps(LHS_matrix, LHS_matrix_init,K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=save_phi, apply_full_relaxation=apply_full_relaxation);
        extractSensorData!(sensors,all_etas,χs,save_phi=save_phi)
        return all_etas, all_phis, sensors
    else
        all_etas = solve_all_timesteps(LHS_matrix, LHS_matrix_init,K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=save_phi, apply_full_relaxation=apply_full_relaxation);
        extractSensorData!(sensors,all_etas,χs,save_phi=save_phi)
        return all_etas, sensors
    end 
end;
