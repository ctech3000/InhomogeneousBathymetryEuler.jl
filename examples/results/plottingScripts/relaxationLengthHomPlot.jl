# skip
using JLD2, CairoMakie
using InhomogeneousBathymetryEuler
set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))

data = load("examples/results/plottingScripts/relaxationLengthHomData.jld2")
etas = data["etas"]
energies = data["energies"]
wave = data["wave"]
transs = data["transs"]
time_vec = data["time_vec"]
χss = data["χss"]
extensions=[0,1,2,4,8]
nExt = length(etas)
xss = [transs[idx_e].x.(χss[idx_e])[end:-1:1] for idx_e = 1:nExt]
idc_χs = [1:length(χss[idx_e]) for idx_e = 1:nExt]
Dt = data["Dts"][1]
χs_phys = χss[1][461:end]
eta_ana = -1/9.81*[transformedAnalyticPotential_dt(χs_phys,zeros(size(χs_phys)),t,0.3,wave,transs[1]) for t in time_vec] 

fig = plotSurfaceOverTime(eta_ana,time_vec,χs_phys,"eta ana")
plotSurfaceOverTime!(fig,etas[1],time_vec,χss[1],"direct")
for idx_e = 2:nExt
    plotSurfaceOverTime!(fig,etas[idx_e],time_vec,χss[idx_e],"L_RX=$(extensions[idx_e])λ")
end
ax = fig.content[1]
xlims!(ax,(-31.87,0))