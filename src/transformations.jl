#= using Ferrite
include("inputParameters.jl")
include("helpfulFunctions.jl") =#

struct σTransform
    x::Function
    z::Function
    χ::Function
    σ::Function
    tBath::Bathymetry
    b_L::Float64
end

function σTransform(domain::AbstractDomain)
    x(χ_) = domain.b_L*χ_ .+ domain.x_L
    z(χ_,σ_) = σ_.*eval_bath(domain.bath,x(χ_))
    χ(x_) = (x_.-domain.x_L)/domain.b_L
    σ(x_,z_) = z_./eval_bath(domain.bath,x_)
    transformedBathymetryPoints = χ.(domain.bath.points)[end:-1:1]
    tBath = Bathymetry(transformedBathymetryPoints,domain.bath.vals[end:-1:1])
    b_L = domain.b_L
    return σTransform(x,z,χ,σ,tBath,b_L)
end

#= function to compute B, B_tilde and D on quadrature points of every cell as well
    as of every facet belonging to the inflow boundary on which a Neumann BC
    is prescribed.
    B and B_tilde are additional factors resulting from σ-transforming
    derivatives. D is the determinant of the jacobian of the σ-
    transform. =#
function compute_B_D(cellvalues::CellValues, facetvalues::FacetValues, dh::DofHandler, domain::AbstractDomain, trans::σTransform)
    b_L = domain.b_L
    nCells = length(dh.grid.cells)
    nFacets_inflow = length(getfacetset(dh.grid,"right"))
    nFacets_surface = length(getfacetset(dh.grid,"bottom"))
    B_domain = Vector{Vector{Float64}}(undef,nCells)
    B_tilde_domain = Vector{Vector{Float64}}(undef,nCells)
    D_domain = Vector{Vector{Float64}}(undef,nCells)
    D_inflow_boundary = Vector{Vector{Float64}}(undef,nFacets_inflow)
    D_surface_boundary = Vector{Vector{Float64}}(undef,nFacets_surface)

    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        node_coords = getcoordinates(cell)
        n_q_points = getnquadpoints(cellvalues)
        q_point_coords = spatial_coordinate.((cellvalues,), collect(1:n_q_points),(node_coords,))
        B_domain[cellid(cell)] = Vector{Float64}(undef,n_q_points)
        B_tilde_domain[cellid(cell)] = Vector{Float64}(undef,n_q_points)
        D_domain[cellid(cell)] = Vector{Float64}(undef,n_q_points)
        for q_point = 1:n_q_points
            σ_point = q_point_coords[q_point][2]
            χ_point = q_point_coords[q_point][1]
            db_point = eval_bath(trans.tBath,χ_point,1)
            b_point = eval_bath(trans.tBath,χ_point,0)
            B_domain[cellid(cell)][q_point] = σ_point*db_point/b_point
            B_tilde_domain[cellid(cell)][q_point] = B_domain[cellid(cell)][q_point]^2 + (b_L/b_point)^2
            D_domain[cellid(cell)][q_point] = abs(b_L*b_point)
        end
    end

    for (facet_idx, facet) in enumerate(FacetIterator(dh, getfacetset(dh.grid,"right")))
        reinit!(facetvalues, facet)
        node_coords = getcoordinates(facet)
        n_q_points = getnquadpoints(facetvalues)
        #q_point_coords = spatial_coordinate.((facetvalues,), collect(1:n_q_points),(node_coords,))
        D_inflow_boundary[facet_idx] = Vector{Float64}(undef,n_q_points)
        for q_point = 1:n_q_points
            D_inflow_boundary[facet_idx][q_point] = abs(b_L)
        end
    end

    for (facet_idx, facet) in enumerate(FacetIterator(dh, getfacetset(dh.grid,"bottom")))
        reinit!(facetvalues, facet)
        node_coords = getcoordinates(facet)
        n_q_points = getnquadpoints(facetvalues)
        #q_point_coords = spatial_coordinate.((facetvalues,), collect(1:n_q_points),(node_coords,))
        D_surface_boundary[facet_idx] = Vector{Float64}(undef,n_q_points)
        for q_point = 1:n_q_points
            D_surface_boundary[facet_idx][q_point] = abs(b_L)
        end
    end

    return B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary
end