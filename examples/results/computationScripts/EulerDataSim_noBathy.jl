#computes data for the comparison between Euler and SWE against measurements, without bathymetry 
using InhomogeneousBathymetryEuler 
using Ferrite, JLD2


wave,_,_ = IrregWave("examples/irregWaveData_noBathy.jld2",inflowDepth=0.3)
lam = 2*pi/computeWavenumber(wave,0.3)
x_L = 0.0
x_D = 15.0
x_R = x_D + 2*lam
x_RXs = [0,-2*lam,-4*lam]
nRX = length(x_RXs)
T = 60
Dx = 0.006
nt = 4801
time_vec = collect(LinRange(0,T,nt))

sensorss = Sensors[]
etas = Vector{Vector{Vector{Float64}}}(undef,nRX)
domains = AbstractDomain[]

for idx_r = eachindex(x_RXs)
    print("run $idx_r/$nRX...\n")
    x_RX = x_RXs[idx_r]
    nχ = round(Integer,(x_R-x_RX)/Dx)
    timeMethod = BackwardDiff()
    outflow = OutflowBC("Dirichlet")
    bathPoints = collect(LinRange(x_L+x_RX,x_R,nχ+1))
    bath = Bathymetry(bathPoints,-0.3*ones(size(bathPoints))) # gauss bathy 
    if idx_r == 1
        domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    else
        domain = RelaxedDampedDomainProperties(x_RX,x_L,x_D,x_R,bath,wave)
    end
    trans = σTransform(domain)
    sensors = Sensors(domain,trans,time_vec)
    grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
    @show Dχ, Dσ, nχ, nσ

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

    all_etas, sensors = solve_all_timesteps!(sensors,LHS_matrix, LHS_matrix_init, K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary);
    etas[idx_r] = all_etas
    push!(sensorss,sensors)
    push!(domains,domain)

end

jldsave("examples/results/plottingScripts/EulerDataSim_noBathyData.jld2";etas,sensorss,domains,time_vec)