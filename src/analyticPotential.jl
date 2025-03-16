#= using Interpolations
include("inputParameters.jl")
include("transformations.jl") =#

#=  contains functions to compute the analytic velocity potential for a simple wave
    physical or transformed coordinates (eg analyticPotential vs transformedAnalyticPotential)
    derivative wrt physical or transformed coordinates (eg ..._dx vs ..._dχ)
    h is the depth of the water as in Newman, hence it's a positive number so h = -bath
    =#
const global GRAV = 9.81

function computeWavenumber(freq::Real,h::Real)
    k = 1
    for i = 1:1000
        k = freq^2/(GRAV*tanh(k*h))
    end
    return k
end

function computeWavenumber(wave::AbstractWave,h::Real)
    return computeWavenumber(getFreq(wave),h)
end

function analyticPotential(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = getAmp(wave)
    if wave.hasFadeIn
        amp *= wave.fadeIn(t)
    end
    freq = getFreq(wave)
    phase = wave.phase
    if isnan(wave.waveNumberAtInflow)
        k = computeWavenumber(freq,h)
    else
        k = wave.waveNumberAtInflow
    end
    return GRAV*amp/freq*cosh(k*(z+h))/cosh(k*h)*sin(k*x - freq*t - phase)
end

function analyticPotential(x::Real, z::Real, t::Real, h::Real, wave::IrregWave)
    val = 0.0
    for compWave in wave.waveList
        val += analyticPotential(x,z,t,h,compWave)
    end
    if wave.hasFadeIn
        val *= wave.fadeIn(t)
    end
    return val
end

function analyticPotential(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::AbstractWave) where T<:Real
    return analyticPotential.(x,z,(t,),(h,),(wave,))
end

function analyticPotential_dx(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = wave.amp
    if wave.hasFadeIn
        amp *= wave.fadeIn(t)
    end
    freq = wave.freq
    phase = wave.phase
    if isnan(wave.waveNumberAtInflow)
        k = computeWavenumber(freq,h)
    else
        k = wave.waveNumberAtInflow
    end
    return k*GRAV*amp/freq*cosh(k*(z+h))/cosh(k*h)*cos(k*x - freq*t - phase)
end

function analyticPotential_dx(x::Real, z::Real, t::Real, h::Real, wave::IrregWave)
    val = 0.0
    for compWave in wave.waveList
        val += analyticPotential_dx(x,z,t,h,compWave)
    end
    if wave.hasFadeIn
        val *= wave.fadeIn(t)
    end
    return val
end

function analyticPotential_dx(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::AbstractWave) where T<:Real
    return analyticPotential_dx.(x,z,(t,),(h,),(wave,))
end

function analyticPotential_dz(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = wave.amp
    if wave.hasFadeIn
        amp *= wave.fadeIn(t)
    end
    freq = wave.freq
    phase = wave.phase
    if isnan(wave.waveNumberAtInflow)
        k = computeWavenumber(freq,h)
    else
        k = wave.waveNumberAtInflow
    end
    return k*GRAV*amp/freq*sinh(k*(z+h))/cosh(k*h)*sin(k*x - freq*t - phase)
end

function analyticPotential_dz(x::Real, z::Real, t::Real, h::Real, wave::IrregWave)
    val = 0.0
    for compWave in wave.waveList
        val += analyticPotential_dz(x,z,t,h,compWave)
    end
    if wave.hasFadeIn
        val *= wave.fadeIn(t)
    end
    return val
end

function analyticPotential_dz(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::AbstractWave) where T<:Real
    return analyticPotential_dz.(x,z,(t,),(h,),(wave,))
end

function analyticPotential_dt(x::Real, z::Real, t::Real, h::Real, wave::SimpleWave)
    amp = wave.amp
    freq = wave.freq
    phase = wave.phase
    if isnan(wave.waveNumberAtInflow)
        k = computeWavenumber(freq,h)
    else
        k = wave.waveNumberAtInflow
    end
    pot_dt = -GRAV*amp*cosh(k*(z+h))/cosh(k*h)*cos(k*x - freq*t - phase)
    if wave.hasFadeIn
        tau = 2*pi/freq
        fadeIn_dt = t < 2*tau ? pi*sin(pi*t/(2*tau))/(4*tau) : 0
        pot = analyticPotential(x,z,t,h,wave)
        return fadeIn_dt*pot + wave.fadeIn(t)*pot_dt
    else
        return pot_dt
    end
end

function analyticPotential_dt(x::Real, z::Real, t::Real, h::Real, wave::IrregWave)
    pot_dt = 0.0
    for compWave in wave.waveList
        pot_dt += analyticPotential_dt(x,z,t,h,compWave)
    end
    if wave.hasFadeIn
        tau = 2*pi/getFreq(wave)
        fadeIn_dt = t < 2*tau ? pi*sin(pi*t/(2*tau))/(4*tau) : 0
        pot = analyticPotential(x,z,t,h,wave)
        return fadeIn_dt*pot + wave.fadeIn(t)*pot_dt
    else
        return pot_dt
    end    
end

function analyticPotential_dt(x::Vector{T}, z::Vector{T}, t::Real, h::Real, wave::AbstractWave) where T<:Real
    return analyticPotential_dt.(x,z,(t,),(h,),(wave,))
end

function transformedAnalyticPotential(χ::Real, σ::Real, t::Real, h::Real, wave::SimpleWave, trans::σTransform)
    x = trans.x(χ)
    z = trans.z(χ,σ)
    return analyticPotential(x,z,t,h,wave)
end

function transformedAnalyticPotential(χ::Real, σ::Real, t::Real, h::Real, wave::IrregWave, trans::σTransform)
    val = 0.0
    for compWave in wave.waveList
        val += transformedAnalyticPotential(χ, σ, t, h, compWave, trans)
    end
    if wave.hasFadeIn
        val *= wave.fadeIn(t)
    end
    return val
end

function transformedAnalyticPotential(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::AbstractWave, trans::σTransform) where T<:Real
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

function transformedAnalyticPotential_dχ(χ::Real, σ::Real, t::Real, h::Real, wave::IrregWave, trans::σTransform)
    val = 0.0
    for compWave in wave.waveList
        val += transformedAnalyticPotential_dχ(χ, σ, t, h, compWave, trans)
    end
    if wave.hasFadeIn
        val *= wave.fadeIn(t)
    end
    return val
end

function transformedAnalyticPotential_dχ(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::AbstractWave, trans::σTransform) where T<:Real
    return transformedAnalyticPotential_dχ.(χ,σ,(t,),(h,),(wave,),(trans,))
end

function transformedAnalyticPotential_dσ(χ::Real, σ::Real, t::Real, h::Real, wave::SimpleWave, trans::σTransform)
    x = trans.x.(χ)
    z = trans.z.(χ,σ)
    dz_part = analyticPotential_dz(x,z,t,h,wave)
    b_part = eval_bath(trans.tBath,χ,1)
    return dz_part*b_part
end

function transformedAnalyticPotential_dσ(χ::Real, σ::Real, t::Real, h::Real, wave::IrregWave, trans::σTransform)
    val = 0.0
    for compWave in wave.waveList
        val += transformedAnalyticPotential_dσ(χ, σ, t, h, compWave, trans)
    end
    if wave.hasFadeIn
        val *= wave.fadeIn(t)
    end
    return val
end

function transformedAnalyticPotential_dσ(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::AbstractWave, trans::σTransform) where T<:Real
    return transformedAnalyticPotential_dσ.(χ,σ,(t,),(h,),(wave,),(trans,))
end

function transformedAnalyticPotential_dt(χ::Real, σ::Real, t::Real, h::Real, wave::SimpleWave, trans::σTransform)
    x = trans.x.(χ)
    z = trans.z.(χ,σ)
    return analyticPotential_dt(x,z,t,h,wave)
end

function transformedAnalyticPotential_dt(χ::Real, σ::Real, t::Real, h::Real, wave::IrregWave, trans::σTransform)
    pot_dt = 0.0
    for compWave in wave.waveList
        pot_dt += transformedAnalyticPotential_dt(χ, σ, t, h, compWave, trans)
    end
    if wave.hasFadeIn
        tau = 2*pi/getFreq(wave)
        fadeIn_dt = t < 2*tau ? pi*sin(pi*t/(2*tau))/(4*tau) : 0
        pot = transformedAnalyticPotential(χ,σ,t,h,wave,trans)
        return fadeIn_dt*pot + wave.fadeIn(t)*pot_dt
    else
        return pot_dt
    end   
end

function transformedAnalyticPotential_dt(χ::Vector{T}, σ::Vector{T}, t::Real, h::Real, wave::AbstractWave, trans::σTransform) where T<:Real
    return transformedAnalyticPotential_dt.(χ,σ,(t,),(h,),(wave,),(trans,))
end