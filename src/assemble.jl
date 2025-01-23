#= using Ferrite
include("inputParameters.jl")
include("domainDiscretization.jl")
include("transformations.jl")
include("analyticPotential.jl") =#

#= function to set dirichlet bc from discretized function on nodes =#
function dirichlet_from_discretized_data(grid::Ferrite.AbstractGrid, field::Symbol, boundary_name::String, dirichlet_data::Vector{Float64})
    boundary = getfacetset(grid, boundary_name)
    boundary_coordinates = get_boundary_coordinates(grid,boundary_name)
    # check, whether boundary is horizontal or vertical
    if abs(boundary_coordinates[2][1] - boundary_coordinates[1][1]) > abs(boundary_coordinates[2][2] - boundary_coordinates[1][2])
        dir = [1 0]
    else
        dir = [0 1]
    end
    interp_points_1D = [dot(point2d,dir) for point2d in boundary_coordinates]
    interpolated_data = linear_interpolation(interp_points_1D, dirichlet_data)
    dbc = Dirichlet(field, boundary, x -> interpolated_data(dot(x,dir)))

    return dbc
end

function assemble_K_element!(Ke::Matrix, cellvalues::CellValues, Be::Vector, B_tilde_e::Vector, De::Vector)
    n_basefuncs = getnbasefunctions(cellvalues)
    fill!(Ke, 0)
    for q_point in 1:getnquadpoints(cellvalues)
        dΩ = getdetJdV(cellvalues, q_point)
        B_point = Be[q_point]
        B_tilde_point = B_tilde_e[q_point]
        D_point = De[q_point]
        for j in 1:n_basefuncs
            ψj_dχ, ψj_dσ  = shape_gradient(cellvalues, q_point, j)
            for i in 1:n_basefuncs
                ψi_dχ, ψi_dσ = shape_gradient(cellvalues, q_point, i)
                Ke[j, i] += (ψi_dχ*ψj_dχ - B_point*ψi_dχ*ψj_dσ - B_point*ψi_dσ*ψj_dχ  + B_tilde_point*ψi_dσ*ψj_dσ)*D_point*dΩ
            end
        end
    end
    return Ke
end

function assemble_K_global(cellvalues::CellValues, dh::DofHandler, domain::AbstractDomain, B_domain::Vector{Vector{Float64}}, B_tilde_domain::Vector{Vector{Float64}}, D_domain::Vector{Vector{Float64}})
    K = allocate_matrix(dh)
    n_basefuncs = getnbasefunctions(cellvalues)
    Ke = zeros(n_basefuncs, n_basefuncs)
    assembler = start_assemble(K)
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        Be = B_domain[cellid(cell)]
        B_tilde_e = B_tilde_domain[cellid(cell)]
        De = D_domain[cellid(cell)]
        assemble_K_element!(Ke, cellvalues, Be, B_tilde_e, De)
        assemble!(assembler, celldofs(cell), Ke)
    end
    K *= 1/domain.b_L^2
    return K
end

function assemble_f_element!(fe::Vector, facetvalues::FacetValues, facet::FacetCache, De::Vector, trans::σTransform, wave::AbstractWave, t_p::Real)
    n_basefuncs = getnbasefunctions(facetvalues)
    fill!(fe,0)
    node_coords = getcoordinates(facet)
    n_q_points = getnquadpoints(facetvalues)
    q_point_coords = spatial_coordinate.((facetvalues,), collect(1:n_q_points),(node_coords,))
    for q_point in 1:n_q_points
        dS = getdetJdV(facetvalues, q_point)
        D_point = De[q_point]
        χ = q_point_coords[q_point][1]
        σ = q_point_coords[q_point][2]
        h = -eval_bath(trans.tBath,χ,0)
        phi_g_dx = analyticPotential_dx(trans.x(χ),trans.z(χ,σ),t_p,h,wave)
        #phi_g_dx = 2

        for j in 1:n_basefuncs
            ψj = shape_value(facetvalues, q_point, j)
            fe[j] += -(ψj*phi_g_dx*D_point)*dS
        end
    end
end

function insert_into_f!(f::Vector, fe::Vector, facetvalues::FacetValues,facet::FacetCache)
    dofs = celldofs(facet)
    for i in 1:getnbasefunctions(facetvalues)
        f[dofs[i]] += fe[i]
    end
end

function assemble_f_global(facetvalues::FacetValues, dh::DofHandler, D_inflow_boundary::Vector{Vector{Float64}}, trans::σTransform, wave::AbstractWave, t_p::Real)
    n_basefuncs = getnbasefunctions(facetvalues)
    fe = zeros(n_basefuncs)
    f = zeros(ndofs(dh))
    for (facet_idx,facet) in enumerate(FacetIterator(dh,getfacetset(dh.grid,"right")))
        reinit!(facetvalues,facet)
        De = D_inflow_boundary[facet_idx]
        assemble_f_element!(fe, facetvalues, facet, De, trans, wave, t_p)
        insert_into_f!(f,fe,facetvalues,facet)
    end
    return f
end

function init_K_f(cellvalues::CellValues, facetvalues::FacetValues, dh::DofHandler,domain::AbstractDomain,B_domain::Vector{Vector{Float64}}, B_tilde_domain::Vector{Vector{Float64}}, D_domain::Vector{Vector{Float64}}, D_inflow_boundary::Vector{Vector{Float64}}, trans::σTransform, wave::AbstractWave, χs::Vector{T}) where T<:Real
    K = assemble_K_global(cellvalues,dh,domain,B_domain,B_tilde_domain,D_domain)
    f = assemble_f_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,0.0)
    K_init = deepcopy(K)

    ch = ConstraintHandler(dh)
    phi_surface_curr1 = transformedAnalyticPotential(χs,0*χs,0.0,-domain.b_L,domain.wave,trans)
    dbc = dirichlet_from_discretized_data(dh.grid, :phi, "bottom", phi_surface_curr1) # "bottom", because in transformed domain coordinates are flipped
    add!(ch, dbc);
    close!(ch)
    apply!(K,f,ch)

    return K, K_init, ch
end

#= function to apply Dirichlet data on prescribed dofs. 
    Doesn't alter K, since the corresponding rows/cols are already zeroed
    after first use of "apply!"
    (code stolen from ConstraintHandler.jl) =#
function apply_dirichlet!(f::Vector,K::SparseMatrixCSC,dirichlet_data::Vector,ch::ConstraintHandler)
    m = meandiag(K)
    for i in 1:length(dirichlet_data)
        d = ch.prescribed_dofs[i]
        v = dirichlet_data[i]
        if v != 0
            for j in nzrange(K, d)
                r = K.rowval[j]
                f[r] -= v * K.nzval[j]
            end
        end
    end
    #= for prescribed dof j we have K[j,j] = meandiag(K_init), so f needs 
    to be adjusted =#
    for i in 1:length(dirichlet_data)
        d = ch.prescribed_dofs[i]
        if length(f) != 0
            vz = dirichlet_data[i]
            f[d] = vz * m
        end
    end
end

function meandiag(K::AbstractMatrix)
    z = zero(eltype(K))
    for i in 1:size(K, 1)
        z += abs(K[i, i])
    end
    return z / size(K, 1)
end