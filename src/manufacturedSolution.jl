struct TrueSolution
    u_0::Function
    u_0_dx::Function
    u_0_dz::Function
    f_0::Function
end

function assemble_f_global(cellvalues::CellValues, dh::DofHandler, B_domain::Vector{Vector{Float64}}, B_tilde_domain::Vector{Vector{Float64}}, D_domain::Vector{Vector{Float64}}, trans::σTransform, u::TrueSolution)
    n_basefuncs = getnbasefunctions(cellvalues)
    fe = zeros(n_basefuncs)
    f = zeros(ndofs(dh))
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        Be = B_domain[cellid(cell)]
        B_tilde_e = B_tilde_domain[cellid(cell)]
        De = D_domain[cellid(cell)]
        assemble_f_element!(fe, cellvalues, cell, De, trans, u)
        insert_into_f!(f, fe, cellvalues, cell)
    end
    return f
end

function assemble_f_element!(fe::Vector, cellvalues::CellValues, cell::CellCache, De::Vector, trans::σTransform, u::TrueSolution)
    n_basefuncs = getnbasefunctions(cellvalues)
    fill!(fe, 0)
    n_q_points = getnquadpoints(cellvalues)
    node_coords = getcoordinates(cell)

    for q_point in 1:n_q_points
        dΩ = getdetJdV(cellvalues, q_point)
        D_point = De[q_point]
        q_point_coords = spatial_coordinate(cellvalues, q_point, node_coords)
        x = trans.x(q_point_coords[1])
        z = trans.z(q_point_coords...)
        f0_point = u.f_0.(x, z)
        for i in 1:n_basefuncs
            Ni = shape_value(cellvalues, q_point, i)
            fe[i] += Ni * f0_point * D_point * dΩ
        end
    end
end

function insert_into_f!(f::Vector, fe::Vector, cellvalues::CellValues, cell::CellCache)
    dofs = celldofs(cell)
    for i in 1:getnbasefunctions(cellvalues)
        f[dofs[i]] += fe[i]
    end
end

function assemble_g_global(facetvalues::FacetValues, dh::DofHandler, D_inflow_boundary::Vector{Vector{Float64}}, trans::σTransform, u::TrueSolution)
    n_basefuncs = getnbasefunctions(facetvalues)
    ge = zeros(n_basefuncs)
    gg = zeros(ndofs(dh))
    for (facet_idx, facet) in enumerate(FacetIterator(dh, getfacetset(dh.grid, "right")))
        reinit!(facetvalues, facet)
        De = D_inflow_boundary[facet_idx]
        assemble_g_element!(ge, facetvalues, facet, De, trans, u)
        insert_into_f!(gg, ge, facetvalues, facet)
    end
    return gg
end

function assemble_g_element!(ge::Vector, facetvalues::FacetValues, facet::FacetCache, De::Vector, trans::σTransform, u::TrueSolution)
    n_basefuncs = getnbasefunctions(facetvalues)
    fill!(ge, 0)
    node_coords = getcoordinates(facet)
    n_q_points = getnquadpoints(facetvalues)

    for q_point in 1:n_q_points
        dS = getdetJdV(facetvalues, q_point)
        D_point = De[q_point]
        q_point_coords = spatial_coordinate(facetvalues, q_point, node_coords)
        x = trans.x(q_point_coords[1])
        z = trans.z(q_point_coords...)
        g_point = -u.u_0_dx(x, z)
        for i in 1:n_basefuncs
            Ni = shape_value(facetvalues, q_point, i)
            ge[i] += (Ni * g_point * D_point) * dS
        end
    end
end

function assemble_h_global(facetvalues::FacetValues, dh::DofHandler, trans::σTransform, u::TrueSolution)
    n_basefuncs = getnbasefunctions(facetvalues)
    he = zeros(n_basefuncs)
    h = zeros(ndofs(dh))
    for (facet_idx, facet) in enumerate(FacetIterator(dh, getfacetset(dh.grid, "bottom")))
        reinit!(facetvalues, facet)
        assemble_h_element!(he, facetvalues, facet, trans, u)
        insert_into_f!(h, he, facetvalues, facet)
    end
    return h
end

function assemble_h_element!(he::Vector, facetvalues::FacetValues, facet::FacetCache, trans::σTransform, u::TrueSolution)
    n_basefuncs = getnbasefunctions(facetvalues)
    fill!(he, 0)
    node_coords = getcoordinates(facet)
    n_q_points = getnquadpoints(facetvalues)

    for q_point in 1:n_q_points
        dS = getdetJdV(facetvalues, q_point)
        q_point_coords = spatial_coordinate(facetvalues, q_point, node_coords)
        D_point = abs(trans.tBath.vals[end])
        x = trans.x(q_point_coords[1])
        z = trans.z(q_point_coords...)
        h_point = u.u_0_dz(x, z)
        for i in 1:n_basefuncs
            Ni = shape_value(facetvalues, q_point, i)
            he[i] += (Ni * h_point * D_point) * dS
        end
    end
end

function coefficientVector(dh::Ferrite.DofHandler, func::Function, trans::σTransform)
    #= coeff = zeros(Float64, ndofs(dh))
    for (i, node) in enumerate(dh.grid.nodes)
        χ, σ = get_node_coordinate(node)
        x = trans.x(χ)
        z = trans.z(χ, σ)
        coeff[i] = func(x, z)
    end
    return coeff =#

    coeff = zeros(Float64, ndofs(dh))

    for cell in CellIterator(dh)
        dofs = celldofs(cell)
        node_coords = getcoordinates(cell)
        for (i,dof) in enumerate(dofs)
            x = trans.x(node_coords[i][1])
            z = trans.z(node_coords[i]...)
            coeff[dof] = func(x, z)
        end
    end
    return coeff
end