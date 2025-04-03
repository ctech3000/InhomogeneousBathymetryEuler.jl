using InhomogeneousBathymetryEuler
using Ferrite, JLD2

nχ_base = 10
factors = [2^i for i = 0:5]
bathTypes = ["flat","gauss0.25","gauss0.5","gauss0.75"]
nHeats1 = length(factors)
nHeats2 = length(bathTypes)
errorsL2 = zeros(Float64,(nHeats1,nHeats2))
errorsMax = zeros(Float64,(nHeats1,nHeats2))
Dχs = zeros(Float64,(nHeats1,nHeats2))
max_u_0s = zeros(Float64,(nHeats1,nHeats2))

for idx_f = 1:nHeats1
    for idx_b = 1:nHeats2
        print("running: factor $(idx_f)/$(nHeats1), bath $(idx_b)/$(nHeats2)...\n")
        x_L = 0.0
        x_R = 1
        timeMethod = BackwardDiff()
        nχ = nχ_base*factors[idx_f]
        nσ = nχ
        outflow = OutflowBC("Dirichlet")
        bathPoints = collect(LinRange(x_L,x_R,nχ+1))
        if bathTypes[idx_b] == "flat"
            bath = Bathymetry(bathPoints,-ones(Float64,nχ+1))
        elseif bathTypes[idx_b] == "gauss0.25"
            bath = Bathymetry(bathPoints,"TrueGauss",shift=0.5,bHeight=0.25,depth=-1.0,sigma=-0.1)
        elseif bathTypes[idx_b] == "gauss0.5"
            bath = Bathymetry(bathPoints,"TrueGauss",shift=0.5,bHeight=0.5,depth=-1.0,sigma=-0.1)
        elseif bathTypes[idx_b] == "gauss0.75"
            bath = Bathymetry(bathPoints,"TrueGauss",shift=0.5,bHeight=0.75,depth=-1.0,sigma=-0.1)
        end

        wave = SimpleWave()
        domain = DomainProperties(x_L,x_R,bath,wave)
        trans = σTransform(domain)
        grid, χs, σs, Dχ, Dσ = discretizeTransformedDomain(domain, trans, nχ=nχ, nσ=nσ)
        @show Dχ, Dσ, nχ, nσ
        Dχs[idx_f,idx_b] = Dχ
        time_vec = collect(0:0.05:30) 

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


        print("K, M done\n")
        
        b(x) = eval_bath(bath,x)
        b_prime(x) = eval_bath(bath,x,1)
        b_pprime(x) = eval_bath(bath,x,2)
        u_0(x,z) = sin(x-domain.x_R)*cosh(z-domain.b_L)
        b_prime_with_val(x) = b(x), b_prime(x)

        @implement_gradient b b_prime_with_val

        function u_0_dx_(x,z,u_0::Function)
            X = Tensors.Tensor{1,2,Float64}((x,z))
            return Tensors.gradient(x->u_0(x[1],x[2]),X)[1]
        end

        function u_0_dz_(x,z,u_0::Function)
            X = Tensors.Tensor{1,2,Float64}((x,z))
            return Tensors.gradient(x->u_0(x[1],x[2]),X)[2]
        end

        function f_0_(x,z,u_0::Function)
            X = Tensors.Tensor{1,2,Float64}((x,z))
            return Tensors.laplace(x->u_0(x[1],x[2]),X)
        end

        u_0_dx(x,z) = u_0_dx_(x,z,u_0)
        u_0_dz(x,z) = u_0_dz_(x,z,u_0)
        f_0(x,z) = f_0_(x,z,u_0)

        f_0_mat = [f_0(trans.x(χ),trans.z(χ,σ)) for χ in χs, σ in σs]

        u = TrueSolution(u_0,u_0_dx,u_0_dz,f_0)
        u_0_coefficients = coefficientVector(dh,u.u_0,trans)
        u_0_dz_coefficients = coefficientVector(dh,u.u_0_dz,trans)
        u_0_nodes = evaluate_at_grid_nodes(dh,u_0_coefficients,:phi) 
        f = assemble_f_global(cellvalues, dh, B_domain, B_tilde_domain, D_domain, trans, u)
        g = assemble_g_global(facetvalues, dh, D_inflow_boundary, trans, u)
        h = M_T0*(u_0_dz_coefficients + u_0_coefficients)
        l = assemble_l_global(facetvalues,dh,domain,trans,u)

        ch = ConstraintHandler(dh)
        dbc = dirichlet_from_discretized_data(dh.grid, :phi, "left", zeros(Float64, nσ+1)) # "left", because in transformed domain coordinates are flipped
        add!(ch, dbc);
        close!(ch)
        RHS = -f + g + h + l
        LHS = K_init + M_T0
        apply!(LHS,RHS,ch)
        print("assembly done\n")

        u_0_nodes_mat = reshape(u_0_nodes,(nχ+1,nσ+1))
        max_u_0 = maximum(abs.(u_0_coefficients))
        max_u_0s[idx_f,idx_b] = max_u_0

        u_num_vec = LHS\RHS
        u_num_nodes = evaluate_at_grid_nodes(dh,u_num_vec,:phi)
        u_num_nodes_mat = reshape(u_num_nodes,(nχ+1,nσ+1));
        errorsL2[idx_f,idx_b] = computeError(u_num_nodes_mat,u_0_nodes_mat,Dχ,norm="L2")
        errorsMax[idx_f,idx_b] = computeError(u_num_nodes_mat,u_0_nodes_mat,Dχ,norm="max")
    end
end

jldsave("examples/results/plottingScripts/convManSolData.jld2"; errorsL2, errorsMax, Dχs, max_u_0s)