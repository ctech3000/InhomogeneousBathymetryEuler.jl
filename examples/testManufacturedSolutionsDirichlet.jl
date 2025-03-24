using InhomogeneousBathymetryEuler, Ferrite, GLMakie, SparseArrays, LinearAlgebra, Tensors

x_L = 0.0
x_R = 3
fac = 1
nχ = 50*fac
nσ = 50*fac
timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
bathPoints = collect(LinRange(x_L,x_R,nχ+1))
#bath = Bathymetry(bathPoints,-3*ones(Float64,nχ+1))
#bath = Bathymetry(bathPoints,"Gauss",shift=1.5,bHeight=0.3,depth=-3.0)
#bath = Bathymetry(bathPoints,"Ramp",rampStart=0.0,rampEnd=3.0,rampHeightStart=-3.0,rampHeightEnd=-2.0)
bath = Bathymetry(bathPoints,"TrueGauss",shift=1.5,bHeight=1, depth=-3.0, sigma=-0.2)

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
#= flat bath, f=0 =#
#u_0(x,z) = sin(x-domain.x_R)*cosh(z-domain.b_L)

#= flat bath, f != 0 =#
#u_0(x,z) = sin(x-domain.x_R)*cosh(z-domain.b_L)*sinh(x-domain.x_R)

#= inhom bath, f != 0 =#
b(x) = eval_bath(bath,x)
b_prime(x) = eval_bath(bath,x,1)
b_pprime(x) = eval_bath(bath,x,2)
u_0(x,z) = sin(x-domain.x_R)*cosh(z-domain.b_L)
b_prime_with_val(x) = b(x), b_prime(x)

@implement_gradient b b_prime_with_val

function u_0_dx_(x,z,u_0::Function)
    X = Tensors.Tensor{1,2,Float64}((x,z))
    return Tensors.gradient(x->u_0(x[1],x[2]),X)[1]
end

function u_0_dz_(x,z,u_0::Function)
    X = Tensors.Tensor{1,2,Float64}((x,z))
    return Tensors.gradient(x->u_0(x[1],x[2]),X)[2]
end

function f_0_(x,z,u_0::Function)
    X = Tensors.Tensor{1,2,Float64}((x,z))
    return Tensors.laplace(x->u_0(x[1],x[2]),X)
end

function b_dx_(x,b_::Function)
    X = Tensors.Tensor{1,2,Float64}((x,0))
    return Tensors.gradient(x->b_(x[1],x[2]),X)[1]
end

u_0_dx(x,z) = u_0_dx_(x,z,u_0)
u_0_dz(x,z) = u_0_dz_(x,z,u_0)
f_0(x,z) = f_0_(x,z,u_0)
b_prime_tensors(x) = b_dx_(x,(x,d) -> b(x))

f_0_mat = [f_0(trans.x(χ),trans.z(χ,σ)) for χ in χs, σ in σs]

u = TrueSolution(u_0,u_0_dx,u_0_dz,f_0)
f = assemble_f_global(cellvalues, dh, B_domain, B_tilde_domain, D_domain, trans, u)
g = assemble_g_global(facetvalues, dh, D_inflow_boundary, trans, u)
u_0_coefficients = coefficientVector(dh,u.u_0,trans)
u_0_dz_coefficients = coefficientVector(dh,u.u_0_dz,trans)
u_0_nodes = evaluate_at_grid_nodes(dh,u_0_coefficients,:phi) 
h = M_T0*u_0_dz_coefficients
l = assemble_l_global(facetvalues,dh,domain,trans,u)

ch = ConstraintHandler(dh)
dbc = dirichlet_from_discretized_data(dh.grid, :phi, "left", zeros(Float64, nσ+1)) # "left", because in transformed domain coordinates are flipped
add!(ch, dbc);
close!(ch)
RHS = -f + g + h + l
apply!(K_init,RHS,ch)
print("assembly done\n")

u_0_nodes_mat = reshape(u_0_nodes,(nχ+1,nσ+1))

u_num_vec = K_init\RHS
u_num_nodes = evaluate_at_grid_nodes(dh,u_num_vec,:phi)
u_num_nodes_mat = reshape(u_num_nodes,(nχ+1,nσ+1));

errorL2 = computeError(u_num_nodes_mat,u_0_nodes_mat,Dχ,norm="L2")
errorMax = computeError(u_num_nodes_mat,u_0_nodes_mat,Dχ,norm="max")

xs = zeros(length(σs))
zs = trans.z(zeros(length(σs)),σs)
u_num_dx_L, _ = computeDerivativeOnBoundary(u_num_nodes_mat,Dχ,Dσ,zeros(length(σs)),σs,domain,trans,"right","physical")
f1 = Figure()
ax1 = Axis(f1[1,1])
lines!(ax1,zs,-u_num_dx_L)
lines!(ax1,zs,-u_0_dx.(xs,zs))

f2 = Figure()
ax2 = Axis(f2[1,1])
xs = trans.x(χs)
zs = zeros(length(xs))
_, u_num_dz_T = computeDerivativeOnBoundary(u_num_nodes_mat,Dχ,Dσ,χs,zeros(length(χs)),domain,trans,"bottom","physical")
lines!(ax2,xs,u_num_dz_T)
lines!(ax2,xs,u_0_dz.(xs,zs))

f3 = Figure()
ax3 = Axis(f3[1,1])
xs = trans.x(χs)
zs = eval_bath(bath,xs)
normals = [computeBathNormal(x,domain) for x in xs]
u_num_dx_B, u_num_dz_B = computeDerivativeOnBoundary(u_num_nodes_mat,Dχ,Dσ,χs,ones(length(χs)),domain,trans,"top","physical")
lines!(ax3,xs,[normals[i][1]*u_num_dx_B[i] + normals[i][2]*u_num_dz_B[i] for i = eachindex(xs)])
lines!(ax3,xs,[normals[i][1]*u_0_dx.(xs[i],zs[i]) + normals[i][2]*u_0_dz.(xs[i],zs[i]) for i = eachindex(xs)])

#= rhs_num = InhomogeneousBathymetryEuler.laplace(u_num_nodes_mat, Dχ, Dσ, trans, χs, σs) - f_0_mat[2:end-1,2:end-1]
rhs_0 = InhomogeneousBathymetryEuler.laplace(u_0_nodes_mat, Dχ, Dσ, trans, χs, σs) - f_0_mat[2:end-1,2:end-1] =#
rhs_num = InhomogeneousBathymetryEuler.laplace(u_num_nodes_mat, Dχ, Dσ, trans, χs, σs)
rhs_0 = InhomogeneousBathymetryEuler.laplace(u_0_nodes_mat, Dχ, Dσ, trans, χs, σs)
f4 = Figure()
ax4 = Axis3(f4[1,1],xlabel="χ",ylabel="σ")
surface!(ax4,χs,σs,rhs_num)

#= f4_5 = Figure()
ax4_5 = Axis(f4_5[1,1],xlabel="χ",ylabel="Δϕ - rhs")
σ_sl = Slider(f4_5[2, 1], range = eachindex(σs), startvalue = 1)
rhs_slice = lift(σ_sl.value) do j
    rhs_num[:,j]
end
lines!(ax4_5,χs[2:end-1],rhs_slice) =#

f5 = Figure()
ax5 = Axis3(f5[1,1],xlabel="χ",ylabel="σ")
surface!(ax5,χs,σs,rhs_0)

#= f6 = Figure()
ax6 = Axis(f6[1,1],xlabel="χ",ylabel="Δϕ - rhs")
σ_sl = Slider(f6[2, 1], range = 1:length(σs50)-2, startvalue = 1)
rhs_slice = lift(σ_sl.value) do j
    rhs_num[:,2*j-1]
end
rhs_slice50 = lift(σ_sl.value) do j
    rhs_num50[:,j]
end
lines!(ax6,χs[2:end-1],rhs_slice)
lines!(ax6,χs50[2:end-1],rhs_slice50) =#
