#= include("transformations.jl")
include("helpfulFunctions.jl")
include("inputParameters.jl") =#



function computeDerivativeOnBoundary(sol::Matrix{Float64},Dχ::Real,Dσ::Real,χs::Vector,σs::Vector,domain::AbstractDomain,trans::σTransform,boundary::String,derivative::String)
    b_L = domain.b_L
    if boundary == "top"
        if derivative == "transformed"
            dσ = 1/Dσ*(3/2*sol[:,end] - 2*sol[:,end-1] + 1/2*sol[:,end-2])
            dχ = firstDerivative(sol[:,end],Dχ)
            return dχ, dσ
        elseif derivative == "physical"
            dσ = 1/Dσ*(3/2*sol[:,end] - 2*sol[:,end-1] + 1/2*sol[:,end-2])
            dχ = firstDerivative(sol[:,end],Dχ)
            dx = dχ*1/b_L + dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
            dz = dσ./(eval_bath(domain.bath,trans.x(χs),0))
            return dx, dz
        end
    elseif boundary == "left"
        if derivative == "transformed"
            dχ = 1/Dχ*(-3/2*sol[1,:] + 2*sol[2,:] - 1/2*sol[3,:])
            dσ = firstDerivative(sol[1,:],Dσ)
            return dχ, dσ
        elseif derivative == "physical"
            dχ = 1/Dχ*(-3/2*sol[1,:] + 2*sol[2,:] - 1/2*sol[3,:])
            dσ = firstDerivative(sol[1,:],Dσ)
            dx = dχ*1/b_L + dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
            dz = dσ./(eval_bath(domain.bath,trans.x(χs),0))
            return dx, dz
        end
    elseif boundary == "right"
        if derivative == "transformed"
            dχ = 1/Dχ*(3/2*sol[end,:] - 2*sol[end-1,:] + 1/2*sol[end-2,:])
            dσ = firstDerivative(sol[end,:],Dσ)
            return dχ, dσ
        elseif derivative == "physical"
            dχ = 1/Dχ*(3/2*sol[end,:] - 2*sol[end-1,:] + 1/2*sol[end-2,:])
            dσ = firstDerivative(sol[end,:],Dσ)
            dx = dχ*1/b_L + dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
            dz = dσ./(eval_bath(domain.bath,trans.x(χs),0))
            return dx, dz
        end
    elseif boundary == "bottom"
        if derivative == "transformed"
            dσ = 1/Dσ*(-3/2*sol[:,1] + 2*sol[:,2] - 1/2*sol[:,3])
            dχ = firstDerivative(sol[:,1],Dχ)
            return dχ, dσ
        elseif derivative == "physical"
            dσ = 1/Dσ*(-3/2*sol[:,1] + 2*sol[:,2] - 1/2*sol[:,3])
            dχ = firstDerivative(sol[:,1],Dχ)
            dx = dχ*1/b_L + dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
            dz = dσ./(eval_bath(domain.bath,trans.x(χs),0))
            return dx, dz
        end
    end
end

function laplace(sol::Matrix{Float64}, Dχ::Real, Dσ::Real, domain::AbstractDomain)
    b_L = domain.b_L
    nχ, nσ = size(sol)
    lap = zeros(Float64, (nχ-2,nσ-2))
    for i = 2:nχ-1
        for j = 2:nσ-1
            val = sol[i-1,j]+sol[i+1,j]+sol[i,j-1]+sol[i,j+1]-4*sol[i,j]
            val *= 1/(Dχ*Dσ)*1/b_L^2
            lap[i-1,j-1] = val
        end
    end
    return lap
end