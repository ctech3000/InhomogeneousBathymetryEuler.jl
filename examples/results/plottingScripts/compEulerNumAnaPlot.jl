using JLD2, CairoMakie
using InhomogeneousBathymetryEuler
set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))
function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I[end:-1:1]
end

save_figs = false

data_2 = load("examples/results/plottingScripts/convEuler_noBathyData.jld2")
etas_2 = data_2["etas"]
wave = data_2["wave"]
Dχs = data_2["Dχs"]
trans = data_2["trans"]
Dxs = -trans.x(Dχs)
time_vec = data_2["time_vec"]
χs = data_2["χs1"]
xs = trans.x.(χs)[end:-1:1]
nFacs = length(etas_2)
equi_inds = 225:250

data_4 = load("examples/results/plottingScripts/convEuler_noBathy_longDampData.jld2")
etas_4 = data_4["etas"]

eta_ana = -1/9.81*[transformedAnalyticPotential_dt(χs,zeros(size(χs)),t,0.3,wave,trans) for t in time_vec]
#max_eta = maximum(maximum([abs.(eta_ana[i]) for i = eachindex(eta_ana)]))
max_eta = 0.005
error_where = 1:length(χs)

errorsL2_t_2 = [[computeError(etas_2[idx_f][idx_t][error_where],eta_ana[idx_t][error_where],Dχs[1],norm="L2") for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsL2_2 = [maximum(errorsL2_t_2[idx_f][equi_inds]) for idx_f = 1:nFacs]
maxPos_2 = [Int64[] for idx_f = 1:nFacs]
errorsMax_t_2 = [[computeError(etas_2[idx_f][idx_t][error_where],eta_ana[idx_t][error_where],Dχs[1],norm="max",maxPos=maxPos_2[idx_f]) for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsMax_2 = [maximum(errorsMax_t_2[idx_f][equi_inds]) for idx_f = 1:nFacs]

errorsL2_t_4 = [[computeError(etas_4[idx_f][idx_t][error_where],eta_ana[idx_t][error_where],Dχs[1],norm="L2") for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsL2_4 = [maximum(errorsL2_t_4[idx_f][equi_inds]) for idx_f = 1:nFacs]
maxPos_4 = [Int64[] for idx_f = 1:nFacs]
errorsMax_t_4 = [[computeError(etas_4[idx_f][idx_t][error_where],eta_ana[idx_t][error_where],Dχs[1],norm="max",maxPos=maxPos_4[idx_f]) for idx_t = eachindex(time_vec)] for idx_f = 1:nFacs]
errorsMax_4 = [maximum(errorsMax_t_4[idx_f][equi_inds]) for idx_f = 1:nFacs]

# error in space over time
CairoMakie.activate!()
fig1_ = Figure(size=(700,300))
#ax11_ = Axis(fig1_[1,1], xlabel="t (s)", ylabel="err", yscale=log10, title="L2")
ax12_ = Axis(fig1_[1,1], xlabel="t (s)", ylabel="error", yscale=log10,xticks=WilkinsonTicks(6,k_min=5))
idx_f = 7
t_inds_f1_ = 1:250

#lines!(ax11_,time_vec[t_inds_f1_],errorsL2_t_2[idx_f][t_inds_f1_],label="L_D=2λ")
lines!(ax12_,time_vec[t_inds_f1_],errorsMax_t_2[idx_f][t_inds_f1_],label="L_D=2λ")
#lines!(ax11_,time_vec[t_inds_f1_],errorsL2_t_4[idx_f][t_inds_f1_],label="L_D=4λ")
lines!(ax12_,time_vec[t_inds_f1_],errorsMax_t_4[idx_f][t_inds_f1_],label="L_D=4λ")

#axislegend(ax11_,position=:lt)
axislegend(ax12_,position=:rt)
ylims!(ax12_,(10^-5,10^-2))
if save_figs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\anaErrTime.svg",fig1_)
end
fig1_

# surface over time
using GLMakie
GLMakie.activate!()
fig2_ = Figure()
ax1_ = Axis(fig2_[1,1],xlabel="χ",ylabel="η",title = "η on whole surface")
t_sl = Slider(fig2_[2, 1], range = eachindex(time_vec), startvalue = 1)
eta_num_t_2 = lift(t_sl.value) do t
    etas_2[4][t]
end
eta_num_t_4 = lift(t_sl.value) do t
    etas_4[4][t]
end
eta_ana_t = lift(t_sl.value) do t
    eta_ana[t]
end
lines!(χs,eta_num_t_2,label="L_D=2λ")
lines!(χs,eta_num_t_4,label="L_D=4λ")
lines!(χs,eta_ana_t,label="ana")
ylims!(ax1_,-1.5*0.005,1.5*0.005)
axislegend(ax1_,position=:rb)
fig2_

#error plot
eocsL2_2 = computeEOC(errorsL2_2,Dχs) 
eocsMax_2 = computeEOC(errorsMax_2,Dχs) 
eocsL2_4 = computeEOC(errorsL2_4,Dχs) 
eocsMax_4 = computeEOC(errorsMax_4,Dχs) 

CairoMakie.activate!()
fig1 = Figure(size=(400,230))
ax1 = Axis(fig1[1,1], xlabel="Δχ (m)", ylabel="relative error",xscale=log10,yscale=log10)
lines!(ax1,Dχs,errorsMax_2/max_eta,label=L"L_D=2λ")
lines!(ax1,Dχs,errorsMax_4/max_eta,label=L"L_D=4λ")
axislegend(ax1,position=:lt)
if save_figs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\anaErrPlot.svg",fig1)
end
fig1


#surface plot
#t_inds_f2 = [80,136,193,250]
t_inds_f2 = [51,76,101,125]
fig2 = Figure(size=(700,450))
axs = Axis[]
for i = 1:4
    push!(axs,Axis(fig2[linToCart(i,2,2)...], xlabel="x (m)", ylabel="η (m)",title="t = $(round(time_vec[t_inds_f2[i]],digits=1))s"))
    lines!(axs[i],xs,eta_ana[t_inds_f2[i]][end:-1:1],label="η ana.")
    lines!(axs[i],xs,etas_4[end][t_inds_f2[i]][end:-1:1],label="η num.")
    #axislegend(axs[i],position=:lt)
end
if save_figs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\anaSurfPlot.svg",fig2)
end
fig2

#= fig3 = Figure()
ax3 = Axis(fig3[1,1],xlabel="t (s)", ylabel="η (m)")
lines!(ax3,time_vec,[etas_2[end][t][1] for t = eachindex(etas_2[end])]) =#