using InhomogeneousBathymetryEuler 
using Ferrite, GLMakie

### select domain size ###
#x_RX = 0.0
x_RX = -20.0
x_L = 0.0
x_D = 15.0
x_R = 25.0
time_vec = collect(0:0.025:25)
nχ = round(Integer,(-x_RX+x_R)/0.0125)   # number of cells in χ direction (nσ will be set accordingly)

plotSolution = false

timeMethod = BackwardDiff()     # time integration method 
outflow = OutflowBC("Dirichlet")    # outflow bc, choose "Dirichlet" or "Neumann"
bathPoints = collect(LinRange(x_L,x_R,nχ+1))

### bathymetries ###
bath = Bathymetry(bathPoints,"Flat")                                          # flat bathy
#bath = Bathymetry(bathPoints,"Gauss",shift=10)                                                 # experiment gauss bathy 
#bath = Bathymetry(bathPoints,"Ramp",rampStart=1,rampEnd=5,rampHeightStart=-3,rampHeightEnd=-1) # ramp bathy
#bath = Bathymetry(bathPoints,"TrueGauss",shift=5,bHeight=0.2,depth=-0.3,sigma=-0.1)            # true gauss bathy

### incoming wave ###
#wave = SimpleWave(inflowDepth=0.3)     # plane wave λ=4.78
wave,_,_ = IrregWave("examples/irregWaveData_noBathy.jld2",inflowDepth=0.3) # irreg wave from sensor data
#wave = IrregWave([0.005,0.0005],[2*pi*0.35,10],[0,0.5],inflowDepth=0.3)

### domain ###
#domain = DampedDomainProperties(x_L,x_D,x_R,bath,wave)
domain = RelaxedDampedDomainProperties(x_RX,x_L,x_D,x_R,bath,wave)

### sensor position ###
sensor_pos = [0.0,2.0,4.0,6.0] # experimental sensor position
#sensor_pos = [4.0,8.0,12.0,16.0]

trans = σTransform(domain)
sensors = Sensors(sensor_pos,domain,trans,time_vec,"x") 
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

# compute values for σ-Transform
B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary = compute_B_D(cellvalues,facetvalues,dh,domain,trans)

# assemble FEM matrices
K, K_init, M_T0, M_T1, M_T2, LHS_matrix, LHS_matrix_init, ch = init_K_M(cellvalues, facetvalues, dh,domain,B_domain, B_tilde_domain, D_domain, D_inflow_boundary, D_surface_boundary, trans, outflow, timeMethod, time_vec, nσ);
print("assembly done, compute solution...\n")

# compute solution
all_etas, sensors = solve_all_timesteps!(sensors,LHS_matrix, LHS_matrix_init, K_init, M_T0, M_T1, domain, trans, χs, σs, time_vec, timeMethod, facetvalues, dh, ch, outflow, D_inflow_boundary, save_phi=false);

# plot solution
if plotSolution
    bathy = true  # computed with bathymetry?
    numTimeInds = collect(1:length(time_vec))  # plot over which time?
    fig = plotSensorData(sensors,time_vec,numLabel="Euler",bathy=bathy,numTimeInds=numTimeInds)
    Legend(fig[3,1:2],fig.content[1],orientation=:horizontal)
    fig
end
