#= using Interpolations
include("helpfulFunctions.jl") =#

abstract type AbstractWave end

struct SimpleWave <: AbstractWave
    amp::Float64
    freq::Float64
    phase::Float64
    fadeIn::Function
    hasFadeIn::Bool
end

function SimpleWave() 
    amp = 0.0053740783
    freq = 2.199114857512855
    phase = 2.7
    T = 2*pi/freq
    fadeInFunction(t) = t < 2*T ? 1/2*(1-cos(pi*t/(2*T))) : 1
    return SimpleWave(amp, freq, phase, fadeInFunction, true)
end

function SimpleWave(str::String)
    if str == "noFadeIn"
        return SimpleWave(0.0053740783,2.199114857512855,2.7,x->1*x,false)
    elseif str == "fadeIn"
        return SimpleWave()
    end
end

struct Bathymetry
    points::Vector{Float64}
    vals::Vector{Float64}
    derivative::Vector{Float64}
end

function Bathymetry(points::Vector{Float64},vals::Vector{Float64})
    derivative = firstDerivative(vals,points[2]-points[1])
    return Bathymetry(points,vals,derivative)
end

function Bathymetry(points::Vector{Float64})
    vals = 1/5*exp.(-(points.-5).^2).-1/3
    return Bathymetry(points,vals)
end

function eval_bath(bath::Bathymetry, x::Real, der::Int=0)
    if der == 0
        interp = linear_interpolation(bath.points, bath.vals)
    elseif der == 1
        interp = linear_interpolation(bath.points, bath.derivative)
    end
    return interp(x)
end

function eval_bath(bath::Bathymetry, x::Vector{T}, der::Int=0) where T<:Real
    if der == 0
        interp = linear_interpolation(bath.points, bath.vals)
    elseif der == 1
        interp = linear_interpolation(bath.points, bath.derivative)
    end
    return interp(x)
end

abstract type AbstractDomain end

struct DomainProperties <: AbstractDomain
    x_L::Float64
    x_R::Float64
    bath::Bathymetry
    b_L::Float64
    wave::AbstractWave
end

function DomainProperties(x_L::Float64, x_R::Float64, bath::Bathymetry, wave::AbstractWave)
    b_L = eval_bath(bath,x_L)
    return DomainProperties(x_L,x_R,bath,b_L,wave)
end

# Standarddomain
function DomainProperties(bathymetryPoints::Vector{Float64}, wave::AbstractWave) 
    bath = Bathymetry(bathymetryPoints)
    x_L = 0.0
    x_R = 20.0
    return DomainProperties(x_L,x_R,bath,wave)
end

function DomainProperties(bathymetryPoints::Vector{Float64})
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
    μ_D(x) = x<x_D ? 0.0 : Float64(wave.freq*((x-x_D)/(x_R-x_D))^2)
    return DampedDomainProperties(x_L,x_D,x_R,bath,b_L,wave,μ_D)
end

function DampedDomainProperties(wave::AbstractWave, bathymetryPoints::Vector{Float64})
    bath = Bathymetry(bathymetryPoints::Vector{Float64})
    x_L = 0.0
    x_D = 12.0
    x_R = 20.0
    return DampedDomainProperties(x_L,x_D,x_R,bath,wave)
end

function DampedDomainProperties(bathymetryPoints::Vector{Float64}) 
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