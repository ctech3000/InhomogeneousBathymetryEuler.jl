using InhomogeneousBathymetryEuler 
using Ferrite


x_L = 0.0
x_D = 10.0
x_R = 15.0
timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
bathPoints = collect(LinRange(x_L,x_R,801))
#bathVals = -1.0*ones(Float64,801)
bathVals = bathPoints*2.5 .- x_R*2.5
bath = Bathymetry(bathPoints,bathVals)
wave = SimpleWave()
#domain = DomainProperties(0.0,5*8.5,bath,wave)
domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
trans = σTransform(domain)