using InhomogeneousBathymetryEuler 
using Ferrite, JLD2


x_L = 0.0
x_D = 10.0
T = 50
nχs = 800 .+ [0,200,400,800,1200,1600,4400]
x_R_diff = [0,2.5,5,10,15,20,55]
nσ = 24
nt = 800
nHeats = length(nχs)

etas = Vector{Vector{Vector{Float64}}}(undef,nHeats)
baths = Vector{Bathymetry}(undef,nHeats)
domains = Vector{DampedDomainProperties}(undef,nHeats)
transs = Vector{σTransform}(undef,nHeats)
Dχs = Vector{Float64}(undef,nHeats)
Dσs = Vector{Float64}(undef,nHeats)
Dts = Vector{Float64}(undef,nHeats)

for i = eachindex(nχs)
    print("run $(i)/$(nHeats)...\n")
    nχ = nχs[i]
    x_R = x_D + x_R_diff[i]
    global nσ
    time_vec = collect(LinRange(0,T,nt))
    timeMethod = BackwardDiff()
    outflow = OutflowBC("Dirichlet")
    bathPoints = collect(LinRange(x_L,x_R,nχ+1))
    bathVals = -0.3*ones(Float64,nχ+1)
    bath = Bathymetry(bathPoints,bathVals)
    wave,_,_ = IrregWave("examples/sensorFreqsWith0.jld2")
    domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    trans = σTransform(domain)
    grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ, nσ=nσ)
    Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)
    Dχs[i] = Dχ
    Dσs[i] = Dσ
    Dts[i] = Dt
    baths[i] = bath
    domains[i] = domain
    transs[i] = trans

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
    all_etas = solve_all_timesteps(LHS_matrix, LHS_matrix_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=false);
    #etas[i] = all_etas
    #phis[i] = all_phis
end


wave = SimpleWave()

jldsave("dampingLayerDataIrregEuler.jld2";etas,Dχs,Dσs,Dts,domains,baths,wave,transs)