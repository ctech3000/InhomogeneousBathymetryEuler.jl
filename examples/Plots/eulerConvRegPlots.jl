using InhomogeneousBathymetryEuler, GLMakie, JLD2

inhom = false
true_eta_ana = false

#= if inhom
    filename = "examples/Plots/convergenceDataRegEulerInhom.jld2"
else
    filename = "convergenceDataRegEuler.jld2"
end
d = load(filename)

wave = d["wave"]
Dσs = d["Dσs"]
Dχs = d["Dχs"]
trans = d["trans"]
domain = d["domain"]
etas = d["etas"]
Dts = d["Dts"]
phis = d["phis"]
#orig_phis = phis
time_vec = collect(0:Dts[1]:20)
χs = collect(trans.χ(domain.x_R):Dχs[1]:0)
nDiscs = length(Dσs)
nt = length(time_vec)
Dχ = Dχs[1]
relevant_χs = 100:201
relevant_ts = 90:100
eta_ana = [-1/GRAV*transformedAnalyticPotential_dt(χs,0*χs,t,0.3,wave,trans) for t in time_vec] =#

# compute errors L2 in space
error_phi_L2_time = Vector{Vector{Float64}}(undef,nDiscs-1)
error_eta_L2_time = Vector{Vector{Float64}}(undef,nDiscs-1)
for disc = 1:nDiscs-1
    error_phi_L2_time[disc] = Vector{Float64}(undef,nt)
    error_eta_L2_time[disc] = Vector{Float64}(undef,nt)
    for t = relevant_ts
        error_phi_L2_time[disc][t] = sqrt(Dχ^2*sum((phis[disc][t]-phis[disc+1][t]).^2))
        #error_phi_L2_time[disc][t] = computeError(phis[disc][t],phis[disc+1][t],Dχ,norm="L2")
        if true_eta_ana
            error_eta_L2_time[disc][t] = computeError(etas[disc][t],eta_ana[t],Dχ,norm="L2")
        else
            error_eta_L2_time[disc][t] = computeError(etas[disc][t],etas[disc+1][t],Dχ,norm="L2")
        end
    end
end
error_phi_L2_max = [maximum(error_phi_L2_time[disc]) for disc = 1:nDiscs-1]
error_eta_L2_max = [maximum(error_eta_L2_time[disc]) for disc = 1:nDiscs-1]

# compute errors L2 in space
error_phi_max_time = Vector{Vector{Float64}}(undef,nDiscs-1)
error_eta_max_time = Vector{Vector{Float64}}(undef,nDiscs-1)
for disc = 1:nDiscs-1
    error_phi_max_time[disc] = Vector{Float64}(undef,nt)
    error_eta_max_time[disc] = Vector{Float64}(undef,nt)
    for t = relevant_ts
        error_phi_max_time[disc][t] = computeError(phis[disc][t],phis[disc+1][t],Dχ,norm="max")
        if true_eta_ana
            error_eta_max_time[disc][t] = computeError(etas[disc][t],eta_ana[t],Dχ,norm="max")
        else
            error_eta_max_time[disc][t] = computeError(etas[disc][t],etas[disc+1][t],Dχ,norm="max")
        end
    end
end
error_phi_max_max = [maximum(error_phi_max_time[disc]) for disc = 1:nDiscs-1]
error_eta_max_max = [maximum(error_eta_max_time[disc]) for disc = 1:nDiscs-1]

if inhom
    suffix = " inhom. bathymetry"
else
    suffix = " hom. bathymetry"
end
fig1 = Figure()
ax11 = Axis(fig1[1,1],xlabel="Dχ",ylabel="e_L2",title="L2 Convergence"*suffix,xscale=log10,yscale=log10)
lines!(ax11,Dχs[1:end-1],error_phi_L2_max,label="phi")
lines!(ax11,Dχs[1:end-1],error_eta_L2_max,label="eta")
axislegend(ax11,position=:lt)
ax12 = Axis(fig1[1,2],xlabel="Dχ",ylabel="e_max",title="max Convergence"*suffix,xscale=log10,yscale=log10)
lines!(ax12,Dχs[1:end-1],error_phi_max_max,label="phi")
lines!(ax12,Dχs[1:end-1],error_eta_max_max,label="eta")
axislegend(ax12,position=:lt)
fig1