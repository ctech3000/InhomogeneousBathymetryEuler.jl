using InhomogeneousBathymetryEuler 
using Ferrite


x_L = 0.0
x_D = 10.0
x_R = 20.0
time_vec = collect(0:0.05:50)
nχ = 800
nσ = 12

timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
bathPoints = collect(LinRange(x_L,x_R,nχ+1))
bathVals = -0.3*ones(Float64,nχ+1)
bath = Bathymetry(bathPoints,bathVals)
#wave = SimpleWave()     #λ=4.78
wave, wave_time = IrregWave("examples/sensorFreqs.jld")
domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
trans = σTransform(domain)
sensors = Sensors(domain,trans,time_vec)
grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)

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


K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);

all_etas_irreg, sensors_irreg = solve_all_timesteps!(sensors,LHS_matrix, LHS_matrix_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=false);


wave = SimpleWave()     #λ=4.78
domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
trans = σTransform(domain)
sensors = Sensors(domain,trans,time_vec)
grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)

B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)

K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);

all_etas_reg, sensors_reg = solve_all_timesteps!(sensors,LHS_matrix, LHS_matrix_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=false);
