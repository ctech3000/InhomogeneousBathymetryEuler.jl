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
            dx = dχ*1/b_L - dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
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
            dx = dχ*1/b_L - dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
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
            dx = dχ*1/b_L - dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
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
            dx = dχ*1/b_L - dσ.*trans.z.(χs,σs).*eval_bath(domain.bath,trans.x(χs),1)./(eval_bath(domain.bath,trans.x(χs),0).^2)
            dz = dσ./(eval_bath(domain.bath,trans.x(χs),0))
            return dx, dz
        end
    end
end

function laplace(sol::Matrix{Float64}, Dχ::Real, Dσ::Real, trans::σTransform, χs::Vector{Float64}, σs::Vector{Float64})
    b_L = trans.b_L
    nχ, nσ = size(sol)
    lap = zeros(Float64, (nχ-2,nσ-2))
    dχχ = 1/Dχ^2*[sol[i+1,j]-2*sol[i,j]+sol[i-1,j] for i=2:nχ-1, j=2:nσ-1]
    dσσ = 1/Dσ^2*[sol[i,j+1]-2*sol[i,j]+sol[i,j-1] for i=2:nχ-1, j=2:nσ-1]
    dχσ = 1/(4*Dσ*Dχ)*[sol[i+1,j+1]+sol[i-1,j-1]-sol[i-1,j+1]-sol[i+1,j-1] for i=2:nχ-1,j=2:nσ-1]
    dσ = 1/(2*Dσ)*[sol[i,j+1]-sol[i,j-1] for i=2:nχ-1,j=2:nσ-1]
    for i = 1:nχ-2
        for j = 1:nσ-2
            χ = χs[i+1]
            σ = σs[j+1]
            b = eval_bath(trans.tBath,χ,0)
            db = eval_bath(trans.tBath,χ,1)
            d2b = eval_bath(trans.tBath,χ,2)
            lap[i,j] = 1/b_L^2*dχχ[i,j] - (σ*b*d2b - 2*σ*db^2)/(b_L^2*b^2)*dσ[i,j] - 2*(σ*db)/(b_L^2*b)*dχσ[i,j] + (((σ*db)/(b_L*b))^2+1/b^2)*dσσ[i,j]
        end
    end
    return lap
end