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

#setup dofs
ip = Lagrange{RefQuadrilateral, 1}()
qr = QuadratureRule{RefQuadrilateral}(2)
qr_facet = FacetQuadratureRule{RefQuadrilateral}(2)
cellvalues = CellValues(qr, ip);
facetvalues = FacetValues(qr_facet, ip)
dh = DofHandler(grid)
add!(dh, :phi, ip)
close!(dh);

t_p = 0.0
B_domain, B_tilde_domain, D_domain, D_inflow_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)

K, K_init, ch = init_K_f(cellvalues,facetvalues,dh,domain,B_domain,B_tilde_domain,D_domain,D_inflow_boundary,trans,wave,χs)

time_vec = collect(0:0.05:10)