using InhomogeneousBathymetryEuler 


bathPoints = collect(LinRange(0,10,200))
bathVals = -0.3*ones(Float64,200)
bath = Bathymetry(bathPoints,bathVals)
wave = SimpleWave()
domain = DomainProperties(0.0,10.0,bath,wave)
trans = σTransform(domain)
nχ = 3
nσ = 3
grid, χs, σs, Dχ, Dσ = discretizeTransformedDomain(domain, trans, nχ=nχ, nσ=nσ)