function computeError(u_ana::Union{Matrix{Float64},Vector{Float64}}, u_num::Union{Matrix{Float64},Vector{Float64}}, Dχ::Real;norm::String="L2")
    if norm == "L2"
        return computeErrorL2(u_ana,u_num, Dχ)
    elseif norm == "max"
        return computeErrorMax(u_ana,u_num)
    else
        print("Error in computeError: Invalid norm!\n")
    end
end

function computeError(u_num::Union{Vector{Matrix{Float64}},Vector{Vector{Float64}}},pairs::Vector{Tuple{Int64,Int64}},Dχs::Vector{Float64};kwargs...)
    errors = zeros(Float64,length(pairs))
    for (idx,pair) in enumerate(pairs)
        Dχ = Dχs[pair[1]]
        u1 = u_num[pair[1]]
        u2 = u_num[pair[2]]
        errors[idx] = computeError(u1,u2,Dχ;kwargs...)
    end
    return errors
end

function computeErrorL2(u_ana::Matrix{Float64}, u_num::Matrix{Float64}, Dχ::Real)
    return sqrt.(sum((u_ana - u_num).^2)*Dχ^2)
end

function computeErrorMax(u_ana::Matrix{Float64}, u_num::Matrix{Float64})
    return maximum(abs.(u_ana - u_num))
end

function computeErrorL2(u_ana::Vector{Float64}, u_num::Vector{Float64}, Dχ::Real)
    return sqrt.(sum((u_ana - u_num).^2)*Dχ)
end

function computeErrorMax(u_ana::Vector{Float64}, u_num::Vector{Float64})
    return maximum(abs.(u_ana - u_num))
end

function computeEOC(errors::Vector{Float64},dxs::Vector{Float64})
    N = length(errors) - 1
    eoc = zeros(N)
    for i = 1:N
        eoc[i] = log(errors[i]/errors[i+1])/log(dxs[i]/dxs[i+1])
    end
    return eoc
end