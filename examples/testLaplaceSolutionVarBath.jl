#some tests on whether the program is able to reproduce the almost correct analytic solution for given dirichlet boundary data 
# and an inhomogeneous bathymetry
# result: yes

using Ferrite, SparseArrays, Interpolations
using GLMakie
using InhomogeneousBathymetryEuler

bathPoints = collect(LinRange(0,10,200))
bath = Bathymetry(bathPoints)
wave = SimpleWave("noFadeIn")
domain = DomainProperties(0.0,10.0,bath,wave)
trans = σTransform(domain)
nχ = 200
nσ = 40
grid, χs, σs, Dχ, Dσ = discretizeTransformedDomain(domain, trans, nχ=nχ, nσ=nσ)

#setup dofs
ip = Lagrange{RefQuadrilateral, 1}()
qr = QuadratureRule{RefQuadrilateral}(2)
qr_facet = FacetQuadratureRule{RefQuadrilateral}(4)
cellvalues = CellValues(qr, ip);
facetvalues = FacetValues(qr_facet, ip)
dh = DofHandler(grid)
add!(dh, :phi, ip)
close!(dh);

t_p = 0.0
B_domain, B_tilde_domain, D_domain, D_inflow_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)
K = assemble_K_global(cellvalues,dh,domain,B_domain,B_tilde_domain,D_domain)
g = assemble_g_global(facetvalues,dh,D_inflow_boundary,trans,domain.wave,t_p)

# set up dirichlet boundary conditions
ch = ConstraintHandler(dh)
phi_surface_curr = zeros(size(χs))
for (i,χ) in enumerate(χs)
    phi_surface_curr[i] = transformedAnalyticPotential(χ,0.0,t_p,-eval(trans.tBath,χ),wave,trans)
end
dbc = dirichlet_from_discretized_data(grid, :phi, "bottom", phi_surface_curr) # "bottom", because in transformed domain coordinates are flipped
add!(ch, dbc);
close!(ch)
apply!(K,f,ch)

phi_vec = evaluate_at_grid_nodes(dh,K\f,:phi)
phi_mat = reshape(phi_vec,(nχ+1,nσ+1))


# test solution 


# test whether Dirichlet cond on surface is fulfilled
sol_top = phi_mat[:,1]
ana = transformedAnalyticPotential(χs,0*χs,t_p,2.0,wave,trans)
diff1 = abs.(sol_top - ana)
f_top = Figure()
ax_top = Axis(f_top[1, 1])
lines!(ax_top,χs,diff1)
# checks out

# test whether Neumann = 0 is fulfilled on outflow boundary
dx_outflow, _ = computeDerivativeOnBoundary(phi_mat, Dχ, Dσ, χs[1]*ones(size(σs)), σs, domain, trans, "left", "physical")
diff2 = abs.(dx_outflow)
f_out = Figure()
ax_out = Axis(f_out[1, 1])
lines!(ax_out,σs,diff2)
# mostly checks out, except near surface where unfitting top Dirichlet

# test whether Neumann = 0 is fulfilled on floor
dx_floor, dz_floor = computeDerivativeOnBoundary(phi_mat, Dχ, Dσ, χs, ones(size(χs)), domain, trans, "top", "physical")
xs = trans.x(χs)
normal = Vector{Vector{Float64}}(undef,length(xs))
for (i,x) in enumerate(xs)
    normal[i] = 1/sqrt(1+eval(bath,x,1)^2)*[eval(bath,x,1),-1]
end
dn_floor = [dx_floor[i]*normal[i][1]+dz_floor[i]*normal[i][2] for i=1:length(xs)]
diff3 = abs.(dn_floor)
f_floor = Figure()
ax_floor = Axis(f_floor[1, 1])
lines!(ax_floor,χs,diff3)

# test whether Neumann = phi_g is fulfilled on inflow
dn_inflow, _ = computeDerivativeOnBoundary(phi_mat, Dχ, Dσ, χs[end]*ones(size(σs)), σs, domain, trans, "right", "physical")
xs = trans.x.(χs[end]*ones(size(σs)))
zs = trans.z(χs[end]*ones(size(σs)), σs)
dn_ana_inflow = analyticPotential_dx(xs, zs, t_p, -domain.b_L, wave)
diff4 = abs.(dn_inflow - dn_ana_inflow)
f = Figure()
ax = Axis(f[1, 1])
lines!(ax, zs, dn_inflow, label="num")
lines!(ax, zs, dn_ana_inflow, label="ana")
axislegend(position = :rb)
f
# checks out

# test with finite differences wheter laplace eq. is fulfilled
lap = laplace(phi_mat, Dχ, Dσ, domain)
f1 = Figure()
ax1 = Axis3(f1[1, 1], xlabel = "χ", ylabel = "σ")
surface!(ax1,χs[2:end-1],σs[2:end-1],lap)
f1
# laplace is fulfilled everywhere except at corner outflow/surface

# plot solution
f2 = Figure()
ax2 = Axis3(f2[1, 1], xlabel = "χ", ylabel = "σ")
surface!(ax2,χs,σs,phi_mat)
f2
# analytic solution
phi_ana = zeros(size(phi_mat))
for i = 1:size(phi_ana,1)
    for j = 1:size(phi_ana,2)
        χ = χs[i]
        σ = σs[j]
        val = transformedAnalyticPotential(χ,σ,t_p,-eval(trans.tBath,χ),wave,trans)
        phi_ana[i,j] = val
    end
end
f3 = Figure()
ax3 = Axis3(f3[1, 1], xlabel = "χ", ylabel = "σ")
surface!(ax3,χs,σs,phi_ana)
f3

#= dχ, dσ = computeDerivativeOnBoundary(phi_mat,Dχ,Dσ,χs[end].+0*σs,σs,domain,trans,"right","transformed")
dx, dz = computeDerivativeOnBoundary(phi_mat,Dχ,Dσ,χs[end].+0*σs,σs,domain,trans,"right","physical")

phi_ana_dχ = transformedAnalyticPotential_dχ(χs[1].+0*σs,σs,t_p,domain.b_L,wave,trans)
#= phi_ana_dx = analyticPotential_dx(0.0*trans.z(χs[end].+0*σs,σs),trans.z(χs[end].+0*σs,σs),0.0,domain.b_L,wave)
phi_ana_dz = analyticPotential_dx(0.0*trans.z(χs[end].+0*σs,σs),trans.z(χs[end].+0*σs,σs),0.0,domain.b_L,wave)
 =#

diff = abs.(dχ .- phi_ana_dχ)

f2 = Figure()
ax2 = Axis(f2[1, 1])
lines!(ax2,σs,dχ,label="num")
lines!(ax2,σs,phi_ana_dχ,label="ana")
axislegend(position = :rb)
f2 =#