using InhomogeneousBathymetryEuler 
using Ferrite


x_L = 0.0
x_D = 10.0
timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
bathPoints = collect(LinRange(0,5*8.5,801))
#bathVals = -1.0*ones(Float64,801)
#bath = Bathymetry(bathPoints,bathVals)
bath = Bathymetry(bathPoints)
wave = SimpleWave()
#domain = DomainProperties(0.0,5*8.5,bath,wave)
domain = DampedDomainProperties(0.0,2.5*8.5,5*8.5,bath,wave)
trans = σTransform(domain)
nχ = 800
nσ = 32
grid, χs, σs, Dχ, Dσ = discretizeTransformedDomain(domain, trans, nχ=nχ, nσ=nσ)
time_vec = collect(0:0.05:30)

#setup dofs
ip = Lagrange{RefQuadrilateral, 1}()
qr = QuadratureRule{RefQuadrilateral}(2)
qr_facet = FacetQuadratureRule{RefQuadrilateral}(4)
cellvalues = CellValues(qr, ip);
facetvalues = FacetValues(qr_facet, ip)
dh = DofHandler(grid)
add!(dh, :phi, ip)
close!(dh);

B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)


K, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);

all_etas, all_phis = solve_all_timesteps(LHS_matrix, LHS_matrix_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=true);
