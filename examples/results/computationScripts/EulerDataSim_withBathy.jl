#computes data for the comparison between Euler and SWE against measurements, with bathymetry 
using InhomogeneousBathymetryEuler 
using Ferrite, JLD2


wave_nB,_,_ = IrregWave("examples/irregWaveData_noBathy.jld2",inflowDepth=0.3)
wave_B,_,_ = IrregWave("examples/irregWaveData_withBathy.jld2",inflowDepth=0.3)
waves = [wave_nB,wave_B]
lam = 2*pi/computeWavenumber(wave_nB,0.3)
x_L = 0.0
x_D = 15.0
x_R = x_D + 2*lam
x_RXs = [0,-2*lam,-4*lam]
nRX = length(x_RXs)
nW = length(waves)
T = 60
Dx = 0.006
nt = 4801
time_vec = collect(LinRange(0,T,nt))

sensorss = Matrix{Sensors}(undef,nW,nRX)
etas = Matrix{Vector{Vector{Float64}}}(undef,nW,nRX)
domains = Matrix{AbstractDomain}(undef,nW,nRX)

for idx_w = eachindex(waves)
    if idx_w == 2
        idx_r_vec = [2]
    else
        idx_r_vec = 1:nRX
    end
    wave = waves[idx_w]
    for idx_r = idx_r_vec
        print("wave $idx_w/$(length(waves)), idx_r $idx_r/$nRX...\n")
        x_RX = x_RXs[idx_r]
        nχ = round(Integer,(x_R-x_RX)/Dx)
        timeMethod = BackwardDiff()
        outflow = OutflowBC("Dirichlet")
        bathPoints = collect(LinRange(x_L+x_RX,x_R,nχ+1))
        bath = Bathymetry(bathPoints,"Gauss") # gauss bathy 
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
        etas[idx_w,idx_r] = all_etas
        sensorss[idx_w,idx_r] = sensors
        domains[idx_w,idx_r] = domain

    end
end

jldsave("examples/results/plottingScripts/EulerDataSim_withBathyData.jld2";etas,sensorss,domains,time_vec)
jldsave("examples/results/plottingScripts/EulerDataSim_withBathyData_onlySensors.jld2";sensorss)