using Ferrite, JLD2
using InhomogeneousBathymetryEuler

wave = IrregWave([0.005,0.001,0.0005],[2*pi*0.35,1.3,10],[0,0.3,0.5],inflowDepth=0.3)
lam = 2*pi/computeWavenumber(wave,0.3)

nt = 3850
phys_wavelengths = 2
nχ_per_lam = 230
extensions = [0,1,2,4,8]
nExt = length(extensions)

x_L = 0.0
x_D = phys_wavelengths*lam
x_Rs = x_D.+extensions*lam
T = 50

etas = Vector{Vector{Vector{Float64}}}(undef,nExt)
energies = Vector{Vector{Float64}}(undef,nExt)
Dχs = Vector{Float64}(undef,nExt)
Dσs = Vector{Float64}(undef,nExt)
Dts = Vector{Float64}(undef,nExt)
χss = Vector{Vector{Float64}}(undef,nExt)
domains = Vector{AbstractDomain}(undef,nExt)
transs = Vector{σTransform}(undef,nExt)

for idx_e = eachindex(extensions)
    print("run $(idx_e)/$(nExt)...\n")
    nχ = (phys_wavelengths+extensions[idx_e])*nχ_per_lam
    x_R = x_Rs[idx_e]
    time_vec = collect(LinRange(0,T,nt))

    timeMethod = BackwardDiff()
    outflow = OutflowBC("Dirichlet")
    bathPoints = collect(LinRange(x_L,x_R,nχ+1))
    bath = Bathymetry(bathPoints,-0.3*ones(size(bathPoints))) # gauss bathy 
    if extensions[idx_e] == 0
        domain = DomainProperties(x_L,x_R,bath,wave)
    else
        domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    end
    trans = σTransform(domain)
    sensors = Sensors(domain,trans,time_vec)
    grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
    @show Dχ, Dσ, nχ, nσ
    Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)
    χss[idx_e] = χs
    transs[idx_e] = trans 
    domains[idx_e] = domain

    #setup dofs
    ip = Lagrange{RefQuadrilateral, 1}()
    qr = QuadratureRule{RefQuadrilateral}(2)
    qr_facet = FacetQuadratureRule{RefQuadrilateral}(4)
    cellvalues = CellValues(qr, ip);
    facetvalues = FacetValues(qr_facet, ip)
    dh = DofHandler(grid)
    add!(dh, :phi, ip)
    close!(dh);
    energy=Float64[]

    ### check next
    B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)


    K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);

    all_etas = solve_all_timesteps(LHS_matrix, LHS_matrix_init, K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary,energy=energy);
    etas[idx_e] = all_etas
    energies[idx_e] = energy
end

time_vec = collect(LinRange(0,T,nt))
jldsave("examples/results/plottingScripts/dampingTestShortPhysData.jld2";etas,energies,Dχs,Dσs,Dts,domains,wave,transs,time_vec,χss)