using InhomogeneousBathymetryEuler, Ferrite, GLMakie, SparseArrays, LinearAlgebra, Tensors

nχs = [25,50,100,200,400]
rampDiff = [0.0,0.25,0.5,1.0]
Dχs = zeros(length(nχs))
errL2 = zeros(length(nχs),length(rampDiff))
errMax = zeros(length(nχs),length(rampDiff))
condsK = zeros(length(nχs),length(rampDiff))
condsLHS = zeros(length(nχs),length(rampDiff))

for rampDiff_idx = eachindex(rampDiff)
    for nχ_idx = eachindex(nχs)
        print("rampIdx: $(rampDiff_idx)/$(length(rampDiff)), nχ_idx: $(nχ_idx)/$(length(nχs))\n")
        x_L = 0.0
        x_R = 3
        nχ = nχs[nχ_idx]
        nσ = nχ
        timeMethod = BackwardDiff()
        outflow = OutflowBC("Dirichlet")
        bathPoints = collect(LinRange(x_L,x_R,nχ+1))
        #bath = Bathymetry(bathPoints,-3*ones(Float64,nχ+1))
        #bath = Bathymetry(bathPoints,"Gauss",shift=0.5,bHeight=0.3,depth=-1.0)
        bath = Bathymetry(bathPoints,"Ramp",rampStart=0.0,rampEnd=3.0,rampHeightStart=-3.0,rampHeightEnd=-3.0+rampDiff[rampDiff_idx])
        #bath = Bathymetry(bathPoints,"TrueGauss",shift=1.5,bHeight=1, depth=-3.0)

        wave = SimpleWave()
        domain = DomainProperties(x_L,x_R,bath,wave)
        #domain = DampedDomainProperties(0.0,2.5*8.5,5*8.5,bath,wave)
        trans = σTransform(domain)
        grid, χs, σs, Dχ, Dσ = discretizeTransformedDomain(domain, trans, nχ=nχ, nσ=nσ)
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

        #= inhom bath, f != 0 =#
        b(x) = eval_bath(bath,x)
        b_prime(x) = eval_bath(bath,x,1)
        b_pprime(x) = eval_bath(bath,x,2)
        u_0(x,z) = sin(x-domain.x_R)*cosh(z-b(x))

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
        f = assemble_f_global(cellvalues, dh, B_domain, B_tilde_domain, D_domain, trans, u)
        g = assemble_g_global(facetvalues, dh, D_inflow_boundary, trans, u)
        u_0_coefficients = coefficientVector(dh,u.u_0,trans)
        u_0_dz_coefficients = coefficientVector(dh,u.u_0_dz,trans)
        u_0_nodes = evaluate_at_grid_nodes(dh,u_0_coefficients,:phi) 
        h = M_T0*u_0_dz_coefficients
        #h = assemble_h_global(facetvalues,dh,trans,u)

        ch = ConstraintHandler(dh)
        dbc = dirichlet_from_discretized_data(dh.grid, :phi, "left", zeros(Float64, nσ+1)) # "left", because in transformed domain coordinates are flipped
        add!(ch, dbc);
        close!(ch)
        RHS = -f + g + h
        #apply_dirichlet!(RHS,K_init,zeros(Float64,nσ+1),ch)
        apply!(K_init,RHS,ch)
        print("assembly done\n")

        if nχ <= 50
            print("computing conds...\n")
            condsK[nχ_idx,rampDiff_idx] = cond(collect(K_init))
            condsLHS[nχ_idx,rampDiff_idx] = cond(collect(LHS_matrix))
            print("done\n")
        end
            


        u_0_nodes_mat = reshape(u_0_nodes,(nχ+1,nσ+1))

        u_num_vec = K_init\RHS
        u_num_nodes = evaluate_at_grid_nodes(dh,u_num_vec,:phi)
        u_num_nodes_mat = reshape(u_num_nodes,(nχ+1,nσ+1));
        
        Dχs[nχ_idx] = Dχ
        errL2[nχ_idx,rampDiff_idx] = computeError(u_0_nodes_mat,u_num_nodes_mat,Dχ,norm="L2")
        errMax[nχ_idx,rampDiff_idx] = computeError(u_0_nodes_mat,u_num_nodes_mat,Dχ,norm="max")
    end
end

f1 = Figure()
ax11 = Axis(f1[1,1],xscale=log2,yscale=log2,title="L2 Error",xlabel="Dχ",ylabel="||e||_2")
ax12 = Axis(f1[1,2],xscale=log2,yscale=log2,title="max Error",xlabel="Dχ",ylabel="||e||_∞")
for rampDiff_idx = eachindex(rampDiff)
    lines!(ax11,Dχs,errL2[:,rampDiff_idx],label="rDiff=$(rampDiff[rampDiff_idx])")
    lines!(ax12,Dχs,errMax[:,rampDiff_idx],label="rDiff=$(rampDiff[rampDiff_idx])")
end
axislegend(ax11,position=:lt)
axislegend(ax12,position=:lt)
f1

f2 = Figure()
ax21 = Axis(f2[1,1],title="cond(K)",xlabel="Dχ")
ax22 = Axis(f2[1,2],title="cond(LHS)",xlabel="Dχ")
for rampDiff_idx = eachindex(rampDiff)
    lines!(ax21,Dχs[1:2],condsK[1:2,rampDiff_idx],label="rDiff=$(rampDiff[rampDiff_idx])")
    lines!(ax22,Dχs[1:2],condsLHS[1:2,rampDiff_idx],label="rDiff=$(rampDiff[rampDiff_idx])")
end
axislegend(ax21,position=:rt)
axislegend(ax22,position=:rt)
f2