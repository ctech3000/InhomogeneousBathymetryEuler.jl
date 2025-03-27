using JLD2, CairoMakie
using InhomogeneousBathymetryEuler
set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))

data = load("examples/results/plottingScripts/dampingTestShortPhysData.jld2")
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
eta_ana = -1/9.81*[transformedAnalyticPotential_dt(χss[1],zeros(size(χss[1])),t,0.3,wave,transs[1]) for t in time_vec] 

# surface over time
using GLMakie
GLMakie.activate!()
fig1_ = Figure()
ax1_ = Axis(fig1_[1,1],xlabel="χ",ylabel="η",title = "η on whole surface")
t_sl = Slider(fig1_[2, 1], range = eachindex(time_vec), startvalue = 1)
eta_num_t = [lift(t_sl.value) do t
    etas[idx_e][t][end:-1:1]
end for idx_e = 1:nExt]
eta_ana_t = lift(t_sl.value) do t
    eta_ana[t][end:-1:1]
end
lines!(ax1_,xss[1],eta_ana_t,label="ana")
for idx_e = 3:nExt
    lines!(ax1_,xss[idx_e],eta_num_t[idx_e],label="L_D=$(extensions[idx_e])λ")
end
ylims!(ax1_,-1.5*0.05,1.5*0.05)
axislegend(ax1_,position=:rb)
fig1_

fig1 = Figure()
ax1 = Axis(fig1[1,1],xlabel="t (s)",ylabel="energy (J)")
for idx_e = 1:nExt
    lines!(ax1,time_vec,energies[idx_e],label="L_D=$(extensions[idx_e])λ")
end
axislegend(ax1,position=:lt)
fig1