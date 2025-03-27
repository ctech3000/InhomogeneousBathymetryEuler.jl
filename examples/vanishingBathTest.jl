using InhomogeneousBathymetryEuler 
using Ferrite
using Printf




function runVanishingBathTest(bHeights::Vector{T}) where T<:Real
    total_tests = length(bHeights)
    all_etas_b = Vector{Vector{Vector{Float64}}}(undef,total_tests)
    x_L = 0.0
    x_D = 10.0
    x_R = 20.0
    time_vec = collect(0:0.05:20)
    nχ = 800
    nσ = 12
    global xs = Vector{Float64}(undef,nχ+1)

    timeMethod = BackwardDiff()
    outflow = OutflowBC("Dirichlet")
    bathPoints = collect(LinRange(x_L,x_R,nχ+1))
    for (test_no,bHeight) in enumerate(bHeights)
        print("run $(test_no)/$(total_tests)...\n")
        bath = Bathymetry(bathPoints,"Gauss",shift=5,bHeight=bHeight)
        wave = SimpleWave()     #λ=4.78
        domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
        trans = σTransform(domain)
        grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
        xs = trans.x(χs)[end:-1:1]
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
        
        all_etas = solve_all_timesteps(LHS_matrix, LHS_matrix_init, K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=false);
        all_etas_b[test_no] = [all_etas[t_ind][end:-1:1] for t_ind = eachindex(time_vec)]
    end
    return all_etas_b, time_vec, xs
end
    
#bHeights = [0.0,0.01,0.02,0.03,0.04,0.05,0.1,0.15,0.20]
bHeights = [0.0,0.003,0.006,0.01,0.015,0.02]
all_etas_b, time_vec, xs = runVanishingBathTest(bHeights)

x_ind = 201
all_etas_at_ind = [[all_etas_b[b_ind][t_ind][x_ind] for t_ind = eachindex(time_vec)] for b_ind = eachindex(bHeights)]

fig = Figure()
ax = Axis(fig[1,1],xlabel="t",ylabel=@sprintf "η at x = %.2f" xs[x_ind])
for (b_ind,bHeight) in enumerate(bHeights)
    lines!(ax,time_vec,all_etas_at_ind[b_ind],label="bH=$(bHeight)")
end
axislegend(ax)
fig

GLMakie.activate!()
amp = 0.005
fig2 = Figure()
ax2 = Axis(fig2[1,1],xlabel="χ",ylabel="η",title = "η on whole surface")
time_ind = 1:length(time_vec)
sg = SliderGrid(
    fig2[2,1],
    (label = "Time", range = time_ind, format=t->(@sprintf "%.2f s" time_vec[t]), startvalue=1),
    tellheight = false
)
#t_sl = Slider(fig2[2, 1], range = time_ind, startvalue = 1)
sliderobservables = [s.value for s in sg.sliders]
etas_time = [lift(sliderobservables...) do t...
    all_etas_b[b_ind][t...] 
end for b_ind = eachindex(bHeights)]
for (b_ind,bHeight) in enumerate(bHeights)
    lines!(ax2,xs,etas_time[b_ind],label="bH=$(bHeights[b_ind])")
end
axislegend(ax2)
ylims!(ax2,-1.5*amp,1.5*amp)
rowsize!(fig2.layout,2,Relative(1/5))
fig2

fig3 = Figure()
ax3 = Axis(fig3[1,1],xlabel="bHeight",ylabel="η",title = (@sprintf("max η at t=%.2f",time_vec[t_ind])))
t_sl = Slider(fig3[2, 1], range = time_ind, startvalue = 1)
data = lift(t_sl.value) do t
    [maximum(all_etas_b[b_ind][t]) for b_ind = eachindex(bHeights)]
end
lines!(ax3,bHeights,data)
scatter!(ax3,bHeights,data)
ylims!(ax3,0.9*amp,1.5*amp)
fig3