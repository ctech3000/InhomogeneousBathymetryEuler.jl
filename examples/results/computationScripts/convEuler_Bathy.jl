using Ferrite, JLD2
using InhomogeneousBathymetryEuler

x_L = 0.0
x_D = 10.0
x_R = 20.0
T = 15
wave = IrregWave([0.005,0.001,0.0005],[2*pi*0.35,1.3,10],[0,0.3,0.5],inflowDepth=0.3)

nχ_base = 200
nt_base = 75
facs = [2^i for i = 0:6]
nFacs = length(facs)

etas = Vector{Vector{Vector{Float64}}}(undef,nFacs)
phis = Vector{Vector{Matrix{Float64}}}(undef,nFacs)
Dχs = Vector{Float64}(undef,nFacs)
Dσs = Vector{Float64}(undef,nFacs)
Dts = Vector{Float64}(undef,nFacs)

for i = eachindex(facs)
    print("run $(i)/$(nFacs)...\n")
    nχ = nχ_base*facs[i]
    nt = (nt_base-1)*facs[i] + 1 
    time_vec = collect(LinRange(0,T,nt))

    timeMethod = BackwardDiff()
    outflow = OutflowBC("Dirichlet")
    bathPoints = collect(LinRange(x_L,x_R,nχ+1))
    bath = Bathymetry(bathPoints,"Gauss",shift=2.5) # gauss bathy 
    domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    trans = σTransform(domain)
    sensors = Sensors(domain,trans,time_vec)
    grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
    @show Dχ, Dσ, nχ, nσ
    Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)
    Dχs[i] = Dχ
    Dσs[i] = Dσ
    Dts[i] = Dt

    #setup dofs
    ip = Lagrange{RefQuadrilateral, 1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    qr_facet = FacetQuadratureRule{RefQuadrilateral}(4)
    cellvalues = CellValues(qr, ip);
    facetvalues = FacetValues(qr_facet, ip)
    dh = DofHandler(grid)
    add!(dh, :phi, ip)
    close!(dh);

    ### check next
    B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)


    K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);

    all_etas, all_phis = solve_all_timesteps(LHS_matrix, LHS_matrix_init, K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=(facs[i],facs[i],facs[i]));
    etas[i] = all_etas
    phis[i] = all_phis
end

bathPoints = collect(LinRange(x_L,x_R,nχ_base*facs[end]+1))
bath = Bathymetry(bathPoints,"Gauss",shift=2.5)
domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
trans = σTransform(domain)
time_vec = collect(LinRange(0,T,nt_base))
xs1 = collect(LinRange(x_L,x_D,round(Integer,nχ_base*facs[1]/2+1)))
χs1 = trans.χ.(xs1[end:-1:1])
jldsave("examples/results/plottingScripts/temp_new/convEuler_BathyData.jld2";etas,phis,Dχs,Dσs,Dts,domain,wave,trans,time_vec,xs1,χs1)