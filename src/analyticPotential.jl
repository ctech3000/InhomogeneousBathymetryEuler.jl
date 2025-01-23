#= using Interpolations
include("inputParameters.jl")
include("transformations.jl") =#

#=  contains functions to compute the analytic velocity potential for a simple wave
    physical or transformed coordinates (eg analyticPotential vs transformedAnalyticPotential)
    derivative wrt physical or transformed coordinates (eg ..._dx vs ..._dχ)
    h is the depth of the water as in Newman, hence it's a positive number so h = -bath
    =#
const global g = 9.81

function computeWaveNumber(freq::Real,h::Real)
    k = 1
    for i = 1:1000
        k = freq^2/(g*tanh(k*h))
    end
    return k
end

function analyticPotential(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = wave.amp
    if wave.hasFadeIn
        amp *= wave.fadeIn(t)
    end
    freq = wave.freq
    phase = wave.phase
    k = computeWaveNumber(freq,h)
    return g*amp/freq*cosh(k*(z+h))/cosh(k*h)*sin(k*x - freq*t - phase)
end

function analyticPotential(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::SimpleWave) where T<:Real
    return analyticPotential.(x,z,(t,),(h,),(wave,))
end

function analyticPotential_dx(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = wave.amp
    if wave.hasFadeIn
        amp *= wave.fadeIn(t)
    end
    freq = wave.freq
    phase = wave.phase
    k = computeWaveNumber(freq,h)
    return k*g*amp/freq*cosh(k*(z+h))/cosh(k*h)*cos(k*x - freq*t - phase)
end

function analyticPotential_dx(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::SimpleWave) where T<:Real
    return analyticPotential_dx.(x,z,(t,),(h,),(wave,))
end

function analyticPotential_dz(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = wave.amp
    if wave.hasFadeIn
        amp *= wave.fadeIn(t)
    end
    freq = wave.freq
    phase = wave.phase
    k = computeWaveNumber(freq,h)
    return k*g*amp/freq*sinh(k*(z+h))/cosh(k*h)*sin(k*x - freq*t - phase)
end

function analyticPotential_dz(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::SimpleWave) where T<:Real
    return analyticPotential_dx.(x,z,(t,),(h,),(wave,))
end

function transformedAnalyticPotential(χ::Real, σ::Real, t::Real, h::Real, wave::SimpleWave, trans::σTransform)
    x = trans.x(χ)
    z = trans.z(χ,σ)
    return analyticPotential(x,z,t,h,wave)
end

function transformedAnalyticPotential(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::SimpleWave, trans::σTransform) where T<:Real
    return transformedAnalyticPotential.(χ,σ,(t,),(h,),(wave,),(trans,))
end

function transformedAnalyticPotential_dχ(χ::Real, σ::Real, t::Real, h::Real, wave::SimpleWave, trans::σTransform)
    x = trans.x.(χ)
    z = trans.z.(χ,σ)
    dx_part = analyticPotential_dx(x,z,t,h,wave)
    dz_part = analyticPotential_dz(x,z,t,h,wave)
    db_part = eval_bath(trans.tBath,χ,1)
    return dx_part*trans.b_L + dz_part*σ*db_part
end

function transformedAnalyticPotential_dχ(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::SimpleWave, trans::σTransform) where T<:Real
    return transformedAnalyticPotential_dχ.(χ,σ,(t,),(h,),(wave,),(trans,))
end

function transformedAnalyticPotential_dσ(χ::Real, σ::Real, t::Real, h::Real, wave::SimpleWave, trans::σTransform)
    x = trans.x.(χ)
    z = trans.z.(χ,σ)
    dz_part = analyticPotential_dz(x,z,t,h,wave)
    b_part = eval_bath(trans.tBath,χ,1)
    return dz_part*b_part
end

function transformedAnalyticPotential_dσ(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::SimpleWave, trans::σTransform) where T<:Real
    return transformedAnalyticPotential_dσ.(χ,σ,(t,),(h,),(wave,),(trans,))
end