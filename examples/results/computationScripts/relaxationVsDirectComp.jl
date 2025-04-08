# test amps behind gauss peak at different positions
using Ferrite, JLD2
using InhomogeneousBathymetryEuler

wave = SimpleWave(0.005,2.199114857512855,0,inflowDepth=0.3)
lam = 2*pi/computeWavenumber(wave,0.3)

nt = 4000
phys_wavelengths = 4
damping_wavelengths = 2
relax_wavelengths = [0,1,2,4]
nχ_per_lam = 230
gauss_shifts = collect(0.5:1/30:1)*lam
waveGenMethods = ["direct","rx","rx","rx"]
nShift = length(gauss_shifts)
nWGM = length(relax_wavelengths)

x_L = 0.0
x_D = phys_wavelengths*lam
x_R = x_D + damping_wavelengths*lam
T = 75

etas = Matrix{Vector{Vector{Float64}}}(undef,nShift,nWGM)
energies = Matrix{Vector{Float64}}(undef,nShift,nWGM)
Dχs = Matrix{Float64}(undef,nShift,nWGM)
Dσs = Matrix{Float64}(undef,nShift,nWGM)
Dts = Matrix{Float64}(undef,nShift,nWGM)
χss = Matrix{Vector{Float64}}(undef,nShift,nWGM)
domains = Matrix{AbstractDomain}(undef,nShift,nWGM)
transs = Matrix{σTransform}(undef,nShift,nWGM)

for idx_s = eachindex(gauss_shifts)
    for idx_w = eachindex(waveGenMethods)
        print("run: shift $(idx_s)/$(nShift), method $(idx_w)/$(nWGM)...\n")
        x_RX = x_L-relax_wavelengths[idx_w]*lam
        nχ = (phys_wavelengths+damping_wavelengths+relax_wavelengths[idx_w])*nχ_per_lam
        time_vec = collect(LinRange(0,T,nt))

        timeMethod = BackwardDiff()
        outflow = OutflowBC("Dirichlet")
        bathPoints = collect(LinRange(x_L,x_R,nχ+1))
        bath = Bathymetry(bathPoints,"Gauss",shift=gauss_shifts[idx_s]) # gauss bathy 
        if waveGenMethods[idx_w] == "direct"
            domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
        elseif waveGenMethods[idx_w] == "rx"
            domain = RelaxedDampedDomainProperties(x_RX,x_L,x_D,x_R,bath,wave)
        else
            domain = 0
        end
        trans = σTransform(domain)
        sensors = Sensors(domain,trans,time_vec)
        grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
        @show Dχ, Dσ, nχ, nσ
        Dt = (time_vec[end] - time_vec[1])/(length(time_vec)-1)
        χss[idx_s,idx_w] = χs
        transs[idx_s,idx_w] = trans 
        domains[idx_s,idx_w] = domain

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
        etas[idx_s,idx_w] = all_etas
        energies[idx_s,idx_w] = energy
    end
end

time_vec = collect(LinRange(0,T,nt))
jldsave("examples/results/plottingScripts/relaxationVsDirectCompData.jld2";etas,energies,Dχs,Dσs,Dts,domains,wave,transs,time_vec,χss,gauss_shifts) 

