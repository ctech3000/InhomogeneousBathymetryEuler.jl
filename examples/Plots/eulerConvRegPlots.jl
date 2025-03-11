using InhomogeneousBathymetryEuler

pwd()
d = load("examples/Plots/convergenceDataRegEuler.jld2")

wave = d["wave"]
Dσs = d["Dσs"]
Dχs = d["Dχs"]
trans = d["trans"]
domain = d["domain"]
etas = d["etas"]
Dts = d["Dts"]
phis = d["phis"]
time_vec = collect(0:Dts[1]:20)
χs = collect(trans.χ(domain.x_R):Dχs[1]:0)
nDiscs = length(Dσs)
nt = length(time_vec)
Dχ = Dχs[1]
relevant_χs = 100:201
relevant_ts = 90:100
eta_ana = [-1/GRAV*transformedAnalyticPotential_dt(χs[relevant_χs],0*χs[relevant_χs],t,0.3,wave,trans) for t in time_vec]

# compute errors L2 in space
error_phi_L2_time = Vector{Vector{Float64}}(undef,nDiscs-1)
error_eta_L2_time = Vector{Vector{Float64}}(undef,nDiscs-1)
for disc = 1:nDiscs-1
    error_phi_L2_time[disc] = Vector{Float64}(undef,nt)
    error_eta_L2_time[disc] = Vector{Float64}(undef,nt)
    for t = relevant_ts
        #error_phi_L2_time[disc][t] = computeError(phis[disc][t]relevant_χs,phis[disc+1][t]relevant_χs,Dχ,norm="L2")
        #error_eta_L2_time[disc][t] = computeError(etas[disc][t][relevant_χs],etas[disc+1][t][relevant_χs],Dχ,norm="L2")
        error_eta_L2_time[disc][t] = computeError(etas[disc][t][relevant_χs],eta_ana[t],Dχ,norm="L2")
    end
end
#error_phi_L2_max = [maximum(error_phi_L2_time[disc]) for disc = 1:nDiscs-1]
error_eta_L2_max = [maximum(error_eta_L2_time[disc]) for disc = 1:nDiscs-1]

# compute errors L2 in space
error_phi_max_time = Vector{Vector{Float64}}(undef,nDiscs-1)
error_eta_max_time = Vector{Vector{Float64}}(undef,nDiscs-1)
for disc = 1:nDiscs-1
    error_phi_max_time[disc] = Vector{Float64}(undef,nt)
    error_eta_max_time[disc] = Vector{Float64}(undef,nt)
    for t = relevant_ts
        #error_phi_max_time[disc][t] = computeError(phis[disc][t]relevant_χs,phis[disc+1][t]relevant_χs,Dχ,norm="max")
        #error_eta_max_time[disc][t] = computeError(etas[disc][t][relevant_χs],etas[disc+1][t][relevant_χs],Dχ,norm="max")
        error_eta_max_time[disc][t] = computeError(etas[disc][t][relevant_χs],eta_ana[t],Dχ,norm="max")
    end
end
#error_phi_max_max = [maximum(error_phi_max_time[disc]) for disc = 1:nDiscs-1]
error_eta_max_max = [maximum(error_eta_max_time[disc]) for disc = 1:nDiscs-1]

fig1 = Figure()
ax11 = Axis(fig1[1,1],xlabel="Dχ",ylabel="e_L2",title="L2 Convergence",xscale=log10,yscale=log10)
#lines!(ax11,Dχs[1:end-1],error_phi_L2_max,label="phi")
lines!(ax11,Dχs[1:end-1],error_eta_L2_max,label="eta")
ax12 = Axis(fig1[1,2],xlabel="Dχ",ylabel="e_max",title="max Convergence",xscale=log10,yscale=log10)
#lines!(ax12,Dχs[1:end-1],error_phi_max_max,label="phi")
lines!(ax12,Dχs[1:end-1],error_eta_max_max,label="eta")
fig1