using GLMakie
using InhomogeneousBathymetryEuler
# test solution 
p = 2
t_p = time_vec[p]
phi_mat = all_phis[p]
phi_surface_vec = all_phis[p][:,1]

# plot computed phi_surface 
f_surface = Figure()
ax_surface = Axis(f_surface[1, 1],xlabel="χ",ylabel="ϕ",title="ϕ on surface")
lines!(ax_surface,χs,phi_surface_vec)

# test whether Dirichlet cond on surface is fulfilled
#= sol_top = phi_mat[:,1]
ana = transformedAnalyticPotential(χs,0*χs,t_p,2.0,wave,trans)
diff1 = abs.(sol_top - ana)
f_top = Figure()
ax_top = Axis(f_top[1, 1])
lines!(ax_top,χs,diff1)
# checks out =#

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
    normal[i] = 1/sqrt(1+eval_bath(bath,x,1)^2)*[eval_bath(bath,x,1),-1]
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
ax = Axis(f[1, 1], xlabel = "z", ylabel = "∇ϕ⋅n",title="inflow BC")
lines!(ax, zs, dn_inflow, label="num")
lines!(ax, zs, dn_ana_inflow, label="ana")
axislegend(position = :lt)
f
# checks out

# test with finite differences wheter laplace eq. is fulfilled
lap = InhomogeneousBathymetryEuler.laplace(phi_mat, Dχ, Dσ, domain)
f1 = Figure()
ax1 = Axis3(f1[1, 1], xlabel = "χ", ylabel = "σ", zlabel = "Δϕ",title="Laplace in domain")
surface!(ax1,χs[2:end-1],σs[2:end-1],lap)
f1
# laplace is fulfilled everywhere except at corner outflow/surface

# plot solution
#= f2 = Figure()
ax2 = Axis3(f2[1, 1], xlabel = "χ", ylabel = "σ", zlabel = "ϕ")
surface!(ax2,χs,σs,phi_mat)
f2 =#