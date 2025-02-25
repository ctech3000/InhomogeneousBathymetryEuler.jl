function firstDerivative(v::Vector{T},Δx::Real) where T<:Real
    dvdx_inner = [(v[i+1] - v[i-1])/2 for i = 2:length(v)-1]
    dvdx1 = -1/2*v[3] + 2*v[2] - 3/2*v[1]
    dvdxend = 1/2*v[end-2] - 2*v[end-1] + 3/2*v[end]
    return 1/Δx*vcat(dvdx1,dvdx_inner,dvdxend)
end

function secondDerivative(v::Vector{T},Δx::Real) where T<:Real
    dvdx_inner = [(v[i+1] - 2*v[i] + v[i-1])/2 for i = 2:length(v)-1]
    #= dvdx1 = 2*v[1] - 5*v[2] + 4*v[3] - v[4]
    dvdxend = 2*v[end] - 5*v[end-1] + 4*v[end-2] - v[end-3] =#
    dvdx1 = 1*v[1] - 2*v[2] + 1*v[3]
    dvdxend = 1*v[end] - 2*v[end-1] + 1*v[end-2]
    return 1/Δx^2*vcat(dvdx1,dvdx_inner,dvdxend)
end

function showIfBigDif(val1::Real,val2::Real,tol::Real)
    dif = abs(val1 - val2)
    dif_rel = dif/abs(val1)
    if dif_rel > tol
        @show val1, val2
    end
end

function dofToCoordinate(dh::Ferrite.DofHandler,i::Integer)
    c = 0
    ind = 0
    for cell in CellIterator(dh)
        if i in celldofs(cell)
            c = cellid(cell)
            ind = findall(x->x==i,celldofs(cell))
            return getcoordinates(cell)[ind][1]
        end
    end
end

function q_pointToDof(q_point::Int)
    if q_point == 1
        return 1
    elseif q_point == 2
        return 2
    elseif q_point == 3
        return 4
    elseif q_point == 4
        return 3
    end
end