using InhomogeneousBathymetryEuler, Ferrite, GLMakie, SparseArrays

x_L = 0.0
x_R = 1.0
nχ = 30
nσ = 30
timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
bathPoints = collect(LinRange(x_L,x_R,nχ+1))
bathVals = -1.0*ones(Float64,nχ+1)
bath = Bathymetry(bathPoints,bathVals)
wave = SimpleWave()
domain = DomainProperties(x_L,x_R,bath,wave)
#domain = DampedDomainProperties(0.0,2.5*8.5,5*8.5,bath,wave)
trans = σTransform(domain)
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

u_0(x,z) = cos(x-domain.x_R)*cosh(z-domain.b_L)
u_0_dx(x,z) = -sin(x-domain.x_R)*cosh(z-domain.b_L)
u_0_dz(x,z) = cos(x-domain.x_R)*sinh(z-domain.b_L)
f_0(x,z) = 0.0*x*z

u = TrueSolution(u_0,u_0_dx,u_0_dz,f_0)
f = assemble_f_global(cellvalues, dh, B_domain, B_tilde_domain, D_domain, trans, u)
g = assemble_g_global(facetvalues, dh, D_inflow_boundary, trans, u)
u_0_coefficients = coefficientVector(dh,u.u_0,trans)
u_0_dz_coefficients = coefficientVector(dh,u.u_0_dz,trans)
u_0_nodes = evaluate_at_grid_nodes(dh,u_0_coefficients,:phi) 
h = M_T0*u_0_dz_coefficients

RHS = f + g + h
u_num_vec = K\RHS
u_num_nodes = evaluate_at_grid_nodes(dh,u_num_vec,:phi)

u_0_nodes_mat = reshape(u_0_nodes,(nχ+1,nσ+1))
u_num_nodes_mat = reshape(u_num_nodes,(nχ+1,nσ+1))