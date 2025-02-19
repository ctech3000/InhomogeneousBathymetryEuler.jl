#= using Interpolations
include("helpfulFunctions.jl") =#

abstract type AbstractWave end

struct SimpleWave <: AbstractWave
    amp::Real
    freq::Real
    phase::Real
    fadeIn::Function
    hasFadeIn::Bool
end

function SimpleWave(amp::Real,freq::Real,phase::Real;hasFadeIn=true)
    T = 2*pi/freq

    fadeIn(t) = hasFadeIn ? (t < 2*T ? 1/2*(1-cos(pi*t/(2*T))) : 1.0) : 1.0
    
    return SimpleWave(amp, freq, phase, fadeIn, hasFadeIn)
end

function SimpleWave() 
    amp = 0.0053740783
    freq = 2.199114857512855
    phase = 2.7
    return SimpleWave(amp, freq, phase, hasFadeIn=true)
end

function SimpleWave(str::String)
    if str == "noFadeIn"
        return SimpleWave(0.0053740783,2.199114857512855,2.7,hasFadeIn=false)
    elseif str == "fadeIn"
        return SimpleWave()
    end
end

struct Bathymetry
    points::Vector{Float64}
    vals::Vector{Float64}
    derivative::Vector{Float64}
    vals_func::Interpolations.Extrapolation
    derivative_func::Interpolations.Extrapolation
    second_derivative_func::Interpolations.Extrapolation
end

function Bathymetry(points::Vector{T1},vals::Vector{T2}) where T1 <:Real where T2 <:Real
    vals_func = linear_interpolation(points, vals, extrapolation_bc=Interpolations.Flat())
    derivative = firstDerivative(vals,points[2]-points[1])
    derivative_func = linear_interpolation(points, derivative, extrapolation_bc=Interpolations.Flat())
    second_derivative = firstDerivative(derivative,points[2]-points[1])
    second_derivative_func = linear_interpolation(points, second_derivative, extrapolation_bc=Interpolations.Flat())
    return Bathymetry(points,vals,derivative,vals_func,derivative_func,second_derivative_func)
end

function Bathymetry(points::Vector{T}, type::String; kargs...) where T<:Real
    if type == "Gauss"
        return Bathymetry(points; kargs...)
    else
        print("Error in Bathymetry: type not valid!")
    end
end

function Bathymetry(points::Vector{T}; shift::Real=2.5) where T<:Real
    interp_points = (shift-0.5875) .+ [0.0, 0.0875, 0.1875, 0.2875, 0.3875, 0.4875, 0.5875, 0.6875, 0.7875, 0.8875, 0.9875, 1.0875, 1.175]
    interp_vals = -0.3 .+ [0.0, 0.024, 0.053, 0.0905, 0.133, 0.182, 0.2, 0.182, 0.133, 0.0905, 0.053, 0.024, 0.0]
    bath_interp = BSplineKit.extrapolate(BSplineKit.interpolate(interp_points,interp_vals,BSplineOrder(4)),BSplineKit.Flat())
    vals = bath_interp.(points)
    return Bathymetry(points,vals)
end

function eval_bath(bath::Bathymetry, x::Real, der::Int=0)
    #x=clamp(x,minimum(bath.points),maximum(bath.points))
    if der == 0
        return bath.vals_func(x)
    elseif der == 1
        return bath.derivative_func(x)
    elseif der == 2
        return bath.second_derivative_func(x)
    end
end

function eval_bath(bath::Bathymetry, x::Vector{T}, der::Int=0) where T<:Real
    return eval_bath.((bath,),x,(der,))
end

abstract type AbstractDomain end

struct DomainProperties <: AbstractDomain
    x_L::Float64
    x_R::Float64
    bath::Bathymetry
    b_L::Float64
    wave::AbstractWave
end

function DomainProperties(x_L::Real, x_R::Real, bath::Bathymetry, wave::AbstractWave)
    b_L = eval_bath(bath,x_L)
    return DomainProperties(x_L,x_R,bath,b_L,wave)
end

# Standarddomain
function DomainProperties(bathymetryPoints::Vector{T}, wave::AbstractWave) where T<:Real
    bath = Bathymetry(bathymetryPoints)
    x_L = 0.0
    x_R = 20.0
    return DomainProperties(x_L,x_R,bath,wave)
end

function DomainProperties(bathymetryPoints::Vector{T}) where T<:Real
    wave = SimpleWave()
    return DomainProperties(bathymetryPoints,wave)
end

struct DampedDomainProperties <: AbstractDomain
    x_L::Float64
    x_D::Float64
    x_R::Float64
    bath::Bathymetry
    b_L::Float64
    wave::AbstractWave
    μ_D::Function
end

function DampedDomainProperties(x_L::Float64, x_D::Float64, x_R::Float64, bath::Bathymetry, wave::AbstractWave)
    b_L = eval_bath(bath,x_L)
    L = x_R - x_D
    λ = 2*pi/computeWaveNumber(wave.freq,b_L)
    μ₀ = -GRAV*log(0.5*10^-5)/(2*wave.freq*L)*sqrt(λ/GRAV)
    μ_D(x) = x<x_D ? 0.0 : μ₀*sqrt(GRAV/λ)*(-2*((x-x_D)/L)^3 + 3*((x-x_D)/L)^2)
    #μ_D(x) = x<x_D ? 0.0 : Float64(wave.freq*((x-x_D)/(x_R-x_D))^2)
    return DampedDomainProperties(x_L,x_D,x_R,bath,b_L,wave,μ_D)
end

function DampedDomainProperties(wave::AbstractWave, bathymetryPoints::Vector{Float64})
    bath = Bathymetry(bathymetryPoints)
    x_L = 0.0
    x_D = 12.0
    x_R = 20.0
    return DampedDomainProperties(x_L,x_D,x_R,bath,wave)
end

function DampedDomainProperties(bathymetryPoints::Vector{T}) where T<:Real
    wave = SimpleWave()
    return DampedDomainProperties(wave,bathymetryPoints)
end

abstract type AbstractTimeSteppingMethod end

struct BackwardDiff <: AbstractTimeSteppingMethod end

struct OutflowBC 
    type::String
end

function OutflowBC(val::Int)
    if val == 0
        return OutflowBC("Dirichlet")
    elseif val == 1
        return OutflowBC("Neumann")
    else
        print("Not a valid input for OutflowBC!")
    end
end