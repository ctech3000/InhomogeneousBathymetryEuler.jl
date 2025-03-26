using JLD2, CairoMakie
using InhomogeneousBathymetryEuler

data = load("examples/results/plottingScripts/convEuler_noBathyData.jld2")
phis = data["phis"]
etas = data["etas"]
wave = data["wave"]
Dχs = data["Dχs"]
Dts = data["Dts"]
trans = data["trans"]
time_vec = data["time_vec"]
χs = data["χs1"]
χs = trans.χ.(0:0.1:10)
nFacs = length(phis)
equi_inds = 1:length(time_vec)


eta_ana = -1/9.81*[transformedAnalyticPotential_dt(χs,zeros(size(χs)),t,0.3,wave,trans) for t in time_vec]
errorsL2_t = [[computeError(etas[idx_f][idx_t],eta_ana[idx_t],Dχs[idx_f],norm="L2") for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsL2 = [maximum([errorsL2_t[idx_f][idx_t] for idx_t = equi_inds]) for idx_f = 1:nFacs]
#errorsMax_t = [computeError([phis[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs,norm="max") for idx_t = eachindex(time_vec)]
#errorsMax = [maximum([errorsMax_t[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]

fig1 = Figure()
ax1 = Axis(fig1[1,1], xlabel="t (s)", ylabel="err", title="L2")
for idx_f = 1:nFacs
    #lines(ax1,equi_inds)
end