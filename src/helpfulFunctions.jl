function firstDerivative(v::Vector,Δx)
    dvdx_inner = [(v[i+1] - v[i-1])/2 for i = 2:length(v)-1]
    dvdx1 = -1/2*v[3] + 2*v[2] - 3/2*v[1]
    dvdxend = 1/2*v[end-2] - 2*v[end-1] + 3/2*v[end]
    return 1/Δx*vcat(dvdx1,dvdx_inner,dvdxend)
end

function showIfBigDif(val1::Real,val2::Real,tol::Real)
    dif = abs(val1 - val2)
    dif_rel = dif/abs(val1)
    if dif_rel > tol
        @show val1, val2
    end
end