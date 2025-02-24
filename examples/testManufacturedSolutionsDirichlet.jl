using InhomogeneousBathymetryEuler, Ferrite, GLMakie, SparseArrays, LinearAlgebra

x_L = 0.0
x_R = 1
nχ = 40
nσ = 40
timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
#bathPoints = collect(LinRange(x_L,x_R,nχ+1))
#bathVals = [-4,-3,-2,-3,-4]
#bathVals = vcat(zeros(Float64, round(Integer,nχ/10)).-0.5,-1.5 .+ cos.(bathPoints[1:end-round(Integer,nχ/10)]))
#bathVals = -1*ones(Float64,nχ+1)
bath = Bathymetry(bathPoints,bathVals)
bath = Bathymetry(bathPoints,"Gauss",shift=0.5)
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
K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);


print("K, M done\n")
# flat bath, f=0
#= u_0(x,z) = sin(x-domain.x_R)*cosh(z-domain.b_L)
u_0_dx(x,z) = cos(x-domain.x_R)*cosh(z-domain.b_L)
u_0_dz(x,z) = sin(x-domain.x_R)*sinh(z-domain.b_L)
f_0(x,z) = 0.0*x*z
f_0_mat = [f_0(trans.x(χ),trans.z(χ,σ)) for χ in χs, σ in σs] =#
# flat bath, f != 0
#= u_0(x,z) = sin(x-domain.x_R)*cosh(z-domain.b_L)*sinh(x-domain.x_R)
u_0_dx(x,z) = cos(x-domain.x_R)*cosh(z-domain.b_L)*sinh(x-domain.x_R) + sin(x-domain.x_R)*cosh(z-domain.b_L)*cosh(x-domain.x_R)
u_0_dz(x,z) = sin(x-domain.x_R)*sinh(z-domain.b_L)*sinh(x-domain.x_R)
f_0(x,z) = 2*cos(x-domain.x_R)*cosh(z-domain.b_L)*cosh(x-domain.x_R) + sin(x-domain.x_R)*cosh(z-domain.b_L)*sinh(x-domain.x_R)
f_0_mat = [f_0(trans.x(χ),trans.z(χ,σ)) for χ in χs, σ in σs] =#

# inhom bath, f != 0
b(x) = eval_bath(bath,x)
b_prime(x) = eval_bath(bath,x,1)
b_pprime(x) = eval_bath(bath,x,2)
u_0(x,z) = sin(x-domain.x_R)*cosh(z-b(x))
u_0_dx(x,z) = cos(x-domain.x_R)*cosh(z-b(x)) - sin(x-domain.x_R)*sinh(z-b(x))*b_prime(x)
u_0_dz(x,z) = sin(x-domain.x_R)*sinh(z-b(x))
f_0(x,z) = -sin(x-domain.x_R)*cosh(z-b(x)) - cos(x-domain.x_R)*sinh(z-b(x))*b_prime(x) - cos(x-domain.x_R)*sinh(z-b(x))*b_prime(x) - sin(x-domain.x_R)*sinh(z-b(x))*b_pprime(x) + sin(x-domain.x_R)*cosh(z-b(x))*((b_prime(x))^2) + sin(x-domain.x_R)*cosh(z-b(x))
f_0_mat = [f_0(trans.x(χ),trans.z(χ,σ)) for χ in χs, σ in σs]

#= u_0(x,z) = (x-x_R)^2*(z-b(x))^2
u_0_dx(x,z) = 2*(x-x_R)*(z-b(x))^2 - 2*(x-x_R)^2*(z-b(x))*b_prime(x)
u_0_dz(x,z) = 2*(x-x_R)^2*(z-b(x))
f_0(x,z) = 2*(x-x_R)*(b(x)-z)*((x-x_R)*b_pprime(x)+2*b_prime(x)) + 2*(x-x_R)*b_prime(x)*((x-x_R)*b_prime(x)+b(x)-z) + 2*(b(x)-z)*((x-x_R)*b_prime(x)+b(x)-z) + 2*(x-x_R)^2
 =#
u = TrueSolution(u_0,u_0_dx,u_0_dz,f_0)
f = assemble_f_global(cellvalues, dh, B_domain, B_tilde_domain, D_domain, trans, u)
g = assemble_g_global(facetvalues, dh, D_inflow_boundary, trans, u)
u_0_coefficients = coefficientVector(dh,u.u_0,trans)
u_0_dz_coefficients = coefficientVector(dh,u.u_0_dz,trans)
u_0_nodes = evaluate_at_grid_nodes(dh,u_0_coefficients,:phi) 
h = M_T0*u_0_dz_coefficients
#h = assemble_h_global(facetvalues,dh,trans,u)

ch = ConstraintHandler(dh)
dbc = dirichlet_from_discretized_data(dh.grid, :phi, "left", zeros(Float64, nσ+1)) # "left", because in transformed domain coordinates are flipped
add!(ch, dbc);
close!(ch)
RHS = -f + g + h
apply_dirichlet!(RHS,K_init,zeros(Float64,nσ+1),ch)
print("assembly done\n")

u_0_nodes_mat = reshape(u_0_nodes,(nχ+1,nσ+1))

u_num_vec = K\RHS
u_num_nodes = evaluate_at_grid_nodes(dh,u_num_vec,:phi)
u_num_nodes_mat = reshape(u_num_nodes,(nχ+1,nσ+1));

lap_num = InhomogeneousBathymetryEuler.laplace(u_num_nodes_mat, Dχ, Dσ, domain)
lap_0 = InhomogeneousBathymetryEuler.laplace(u_0_nodes_mat, Dχ, Dσ, domain);