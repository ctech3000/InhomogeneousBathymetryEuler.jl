using JLD2, CairoMakie
using InhomogeneousBathymetryEuler
set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))

data = load("examples/results/plottingScripts/temp_new/convEuler_noBathyData.jld2")
phis = data["phis"]
etas = data["etas"]
wave = data["wave"]
Dχs = data["Dχs"]
Dts = data["Dts"]
trans = data["trans"]
time_vec = data["time_vec"]
χs = data["χs1"]
xs = trans.x.(χs)[end:-1:1]
nFacs = length(phis)
equi_inds = 65:length(time_vec)


eta_ana = -1/9.81*[transformedAnalyticPotential_dt(χs,zeros(size(χs)),t,0.3,wave,trans) for t in time_vec]
errorsL2_t = [[computeError(etas[idx_f][idx_t],eta_ana[idx_t],Dχs[1],norm="L2") for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsL2 = [maximum(errorsL2_t[idx_f][equi_inds]) for idx_f = 1:nFacs]
errorsMax_t = [[computeError(etas[idx_f][idx_t],eta_ana[idx_t],Dχs[1],norm="max") for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsMax = [maximum(errorsMax_t[idx_f][equi_inds]) for idx_f = 1:nFacs]

# error in space over time
CairoMakie.activate!()
fig1_ = Figure()
ax11_ = Axis(fig1_[1,1], xlabel="t (s)", ylabel="err", title="L2")
ax12_ = Axis(fig1_[1,2], xlabel="t (s)", ylabel="err", title="max")
for idx_f = 1:nFacs
    lines!(ax11_,time_vec[equi_inds],errorsL2_t[idx_f][equi_inds],label="2^$(idx_f-1)")
    lines!(ax12_,time_vec[equi_inds],errorsMax_t[idx_f][equi_inds],label="2^$(idx_f-1)")
end
axislegend(ax11_,position=:lt)
fig1_

# surface over time
using GLMakie
GLMakie.activate!()
fig2_ = Figure()
ax1_ = Axis(fig2_[1,1],xlabel="χ",ylabel="η",title = "η on whole surface")
t_sl = Slider(fig2_[2, 1], range = equi_inds, startvalue = 1)
eta_num_t = lift(t_sl.value) do t
    etas[end][t]
end
eta_ana_t = lift(t_sl.value) do t
    eta_ana[t]
end
lines!(χs,eta_num_t,label="num")
lines!(χs,eta_ana_t,label="ana")
ylims!(ax1_,-1.5*0.05,1.5*0.05)
axislegend(ax1_,position=:rb)
fig2_

#error plot
eocsL2 = computeEOC(errorsL2,Dχs) 
eocsMax = computeEOC(errorsMax,Dχs) 

CairoMakie.activate!()
fig1 = Figure(size=(400,230))
ax1 = Axis(fig1[1,1], xlabel="t (s)", ylabel="relative error",xscale=log10,yscale=log10)
lines!(ax1,Dχs,errorsL2,label="L2 error")
lines!(ax1,Dχs,errorsMax,label="Max error")
axislegend(ax1,position=:lt)
fig1

print("\"converged\" grid: Dχ = $(Dχs[5]),  Dt = $(Dts[5])\n")

#surface plot
t_ind = 70
fig2 = Figure(size=(600,450))
ax2 = Axis(fig2[1,1], xlabel="x (m)", ylabel="η",yticks=[-0.05,-0.005,0,0.005,0.05])
lines!(ax2,xs,eta_ana[t_ind][end:-1:1],label="η ana.")
lines!(ax2,xs,etas[end-1][t_ind][end:-1:1],label="η num.")
lines!(ax2,xs,etas[end-1][t_ind][end:-1:1]-eta_ana[t_ind][end:-1:1],label="difference")
axislegend(ax2,position=:lt)
fig2
