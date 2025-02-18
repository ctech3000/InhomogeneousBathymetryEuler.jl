#= using Ferrite
include("inputParameters.jl")
include("transformations.jl")
 =#
function side_idx_to_node_number(side_idx::Int)
    if side_idx == 1
        return [1,2]
    elseif side_idx == 2
        return [2,3]
    elseif side_idx == 3
        return [3,4]
    elseif side_idx == 4
        return [4,1]
    end
end

#= function to compute coordinates of boundary nodes which are
    part of boundary boundary_name =#
function get_boundary_coordinates(grid::Ferrite.AbstractGrid, boundary_name::String)
    n_cells = length(grid.facetsets[boundary_name])
    facet_indices = [face.idx for face in grid.facetsets[boundary_name]]
    coordinates = Vector{Vector{Float64}}(undef,2*n_cells)
    for (idx,face_idx) in enumerate(facet_indices)
        cell_boundary_coordinates = getcoordinates(grid,face_idx[1])[side_idx_to_node_number(face_idx[2])]

        coordinates[2*(idx-1)+1:2*(idx-1)+2] = cell_boundary_coordinates
    end
    unique!(coordinates)
    if boundary_name == "left" || boundary_name == "right"
        sort!(coordinates, by = x -> x[2])
    elseif boundary_name == "top" || boundary_name == "bottom"
        sort!(coordinates, by = x -> x[1])
    else
        print("Error in get_boundary_coordinates: Invalid boundary_name!")
    end
    return coordinates
end

function discretizeTransformedDomain(domain::AbstractDomain, trans::σTransform; nχ::Int=-1, nσ::Int=-1)
    σ_max = 1.0
    σ_min = 0.0
    χ_max = 0.0
    χ_min = trans.χ(domain.x_R)
    if nσ < 0
        if nχ < 0
            print("ERROR: Either nχ or nσ has to be defined!")
        else
            nσ = max(round(Int,abs(σ_max - σ_min)*nχ/abs(χ_max-χ_min)),1)
        end
    elseif nχ < 0
        if nσ < 0
            print("ERROR: Either nχ or nσ has to be defined!")
        else
            nχ = max(round(Int,abs(χ_max - χ_min)*nσ/abs(σ_max-σ_min)),1)
        end
    end
    χ_lims = Ferrite.Vec(χ_min,χ_max)
    σ_lims = Ferrite.Vec(σ_min,σ_max)
    grid = generate_grid(Quadrilateral,(nχ,nσ),χ_lims,σ_lims)
    χs = [point[1] for point in get_boundary_coordinates(grid,"top")]
    σs = [point[2] for point in get_boundary_coordinates(grid,"right")]
    Dχ = χs[2] - χs[1]
    Dσ = σs[2] - σs[1]

    return grid, χs, σs, Dχ, Dσ
end
