using InhomogeneousBathymetryEuler, Ferrite, LinearAlgebra, BenchmarkTools, SparseArrays, IterativeSolvers

function plain(LHS::SparseMatrixCSC{Float64, Int64}, RHS::Vector{Float64})
    return LHS\RHS
end

function withLU(LHS::SparseArrays.UMFPACK.UmfpackLU{Float64, Int64},RHS::Vector{Float64})
    return LHS\RHS
end

function withCG(LHS::SparseMatrixCSC{Float64, Int64}, RHS::Vector{Float64})
    return cg(LHS,RHS)
end

x_L = 0.0
x_D = 10.0
x_R = 20.0
nχ_base = 800
nσ_base = 12
for i = 0:4
    nχ = nχ_base*2^i
    nσ = nσ_base*2^i
    print("nχ = $(nχ)...\n")
    time_vec = collect(0:0.05:20)
    timeMethod = BackwardDiff()
    outflow = OutflowBC("Dirichlet")
    bathPoints = collect(LinRange(x_L,x_R,nχ+1))
    bathVals = -0.3*ones(Float64,nχ+1)
    bath = Bathymetry(bathPoints,bathVals)
    wave = SimpleWave()
    domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
    trans = σTransform(domain)
    sensors = Sensors(domain,trans,time_vec)
    grid, χs, σs, Dχ, Dσ, nχ, nσ = discretizeTransformedDomain(domain, trans, nχ=nχ)
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

    t1 = time()
    K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);
    t2 = time()
    print("Assembly took $(t2-t1)s\n")
    RHS = rand((nχ+1)*(nσ+1))
    if outflow.type == "Dirichlet"
        apply_dirichlet!(RHS,LHS_matrix_init,zeros(Float64,nσ+1),ch)
    end

    print("benchmarking plain...\n")
    @btime plain($LHS_matrix,$RHS)
    t1 = time()
    LU = lu(LHS_matrix)
    t2 = time()
    print("LU factorization took $(t2-t1)s\n")
    print("benchmarking LU...\n")
    @btime withLU($LU,$RHS)
    print("benchmarking CG...\n")
    @btime withCG($LHS_matrix,$RHS)

    x_plain = plain(LHS_matrix,RHS)
    x_LU = withLU(LU,RHS)
    x_CG = withCG(LHS_matrix,RHS)
    print("errors:      LU:$(maximum(abs.(x_plain-x_LU)))       CG:$(maximum(abs.(x_plain-x_CG)))\n")
    print("\n")
end