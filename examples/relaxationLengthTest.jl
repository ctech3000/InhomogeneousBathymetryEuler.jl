using InhomogeneousBathymetryEuler 
using Ferrite, JLD2

wave,amp0,waveTime = IrregWave("examples/irregWaveData_noBathy.jld2",inflowDepth=0.3) # irreg wave
lam = 2*pi/computeWavenumber(wave,0.3)
extension_vec = [0.0,0.5,1.0,1.5,2.0,3.0,4.0]
x_RX_vec = -extension_vec*lam
x_L = 0.0
x_D = 10.0
x_R = 20.0
time_vec = collect(0:0.05/2:3*17.5+30)
nχ_per_λ = 400*2
nχ_vec = round.(Integer,nχ_per_λ*(extension_vec .+ 4))
timeMethod = BackwardDiff()
outflow = OutflowBC("Dirichlet")
nHeats = length(extension_vec)

etas_heats = Vector{Vector{Vector{Float64}}}(undef,nHeats)
sensors_heats = Vector{Sensors}(undef,nHeats)
χs_heats = Vector{Vector{Float64}}(undef,nHeats)
domains_heats = Vector{AbstractDomain}(undef,nHeats)
transs_heats = Vector{σTransform}(undef,nHeats)

for heat = 1:nHeats
    print("heat $(heat)/$(nHeats)...\n")
    nχ = nχ_vec[heat]
    x_RX = x_RX_vec[heat]
    bathPoints = collect(LinRange(x_L,x_R,nχ+1))
    bath = Bathymetry(bathPoints,"Gauss",shift=2.5) # gauss bathy 
    #domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    if x_RX == 0.0
        domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    else
        domain = RelaxedDampedDomainProperties(x_RX,x_L,x_D,x_R,bath,wave)
    end
    trans = σTransform(domain)
    sensors = Sensors(domain,trans,time_vec)
    grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
    Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)
    @show x_RX, Dχ, Dσ, nσ

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

    all_etas, sensors = solve_all_timesteps!(sensors,LHS_matrix, LHS_matrix_init, K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=false);
    etas_heats[heat] = all_etas
    sensors_heats[heat] = sensors
    χs_heats[heat] = χs
    domains_heats[heat] = domain
    transs_heats[heat] = trans
end

jldsave("examples/Plots/relaxationLengthDataFiner.jld2"; etas_heats, sensors_heats, χs_heats, domains_heats, transs_heats, time_vec, extension_vec)