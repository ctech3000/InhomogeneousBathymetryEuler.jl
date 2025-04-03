using JLD2, CairoMakie
using InhomogeneousBathymetryEuler

pairs = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)]
nPairs = length(pairs)
pairsFirst = [pairs[i][1] for i = 1:nPairs]

data_b = load("examples/results/plottingScripts/convEuler_BathyData.jld2")
phis_b = data_b["phis"]
nFacs = length(phis_b)
etas_b = data_b["etas"]
Dχs_b = data_b["Dχs"]
Dts_b = data_b["Dts"]
time_vec = data_b["time_vec"]
errorsL2_t_b = [computeError([phis_b[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_b[1]*ones(size(Dχs_b)),norm="L2") for idx_t = eachindex(time_vec)]
errorsL2_b = [maximum([errorsL2_t_b[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]
errorsMax_t_b = [computeError([phis_b[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_b[1]*ones(size(Dχs_b)),norm="max") for idx_t = eachindex(time_vec)]
errorsMax_b = [maximum([errorsMax_t_b[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]

data_nb = load("examples/results/plottingScripts/convEuler_noBathyData.jld2")
phis_nb = data_nb["phis"]
etas_nb = data_nb["etas"]
Dχs_nb = data_nb["Dχs"]
Dts_nb = data_nb["Dts"]
χs_nb = data_nb["χs1"]
errorsL2_t_nb = [computeError([phis_nb[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_nb[1]*ones(size(Dχs_nb)),norm="L2") for idx_t = eachindex(time_vec)]
errorsL2_nb = [maximum([errorsL2_t_nb[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]
errorsMax_t_nb = [computeError([phis_nb[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_nb[1]*ones(size(Dχs_nb)),norm="max") for idx_t = eachindex(time_vec)]
errorsMax_nb = [maximum([errorsMax_t_nb[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]


set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))
baths = ["Flat bath.","Gauß bath."]
nBaths = length(baths)
errorsL2 = [errorsL2_nb,errorsL2_b]
errorsMax = [errorsMax_nb, errorsMax_b]
Dχs = [Dχs_nb,Dχs_b]

fig1 = Figure(size=(700,350))
ax11 = Axis(fig1[1,1],xlabel="Δχ",ylabel="error",xscale=log10,yscale=log10,title="L2 Error")
for idx_b = 1:nBaths
    lines!(ax11,Dχs_b[pairsFirst],errorsL2[idx_b],label=baths[idx_b])
end
axislegend(ax11,position=:rb)
ax12 = Axis(fig1[1,2],xlabel="Δχ",ylabel="error",xscale=log10,yscale=log10,title="Maximum Error")
for idx_b = 1:nBaths
    lines!(ax12,Dχs_b[pairsFirst],errorsMax[idx_b],label=baths[idx_b])
end
axislegend(ax12,position=:rb)
fig1
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\phiConvErr.svg",fig1)


eocsL2_phi = [computeEOC(errorsL2[idx_b],Dχs[idx_b]) for idx_b = 1:nBaths]
eocsMax_phi = [computeEOC(errorsMax[idx_b],Dχs[idx_b]) for idx_b = 1:nBaths]
print("### PHI ###\nEOCs in L2 and max norm (idx_b x idx_f):\n")
@show eocsL2_phi
@show eocsMax_phi
avgL2 = sum.([eocsL2_phi[idx_b] for idx_b=1:nBaths])/(nFacs-2)
avgMax = sum.([eocsMax_phi[idx_b] for idx_b=1:nBaths])/(nFacs-2)
print("average:\nL2: $(avgL2),  Max: $(avgMax)\n")


### eta

errorsL2_t_b = [computeError([etas_b[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_b[1]*ones(size(Dχs_b)),norm="L2") for idx_t = eachindex(time_vec)]
errorsL2_b = [maximum([errorsL2_t_b[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]
errorsMax_t_b = [computeError([etas_b[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_b[1]*ones(size(Dχs_b)),norm="max") for idx_t = eachindex(time_vec)]
errorsMax_b = [maximum([errorsMax_t_b[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]

errorsL2_t_nb = [computeError([etas_nb[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_nb[1]*ones(size(Dχs_nb)),norm="L2") for idx_t = eachindex(time_vec)]
errorsL2_nb = [maximum([errorsL2_t_nb[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]
errorsMax_t_nb = [computeError([etas_nb[idx_f][idx_t] for idx_f = 1:nFacs],pairs,Dχs_nb[1]*ones(size(Dχs_nb)),norm="max") for idx_t = eachindex(time_vec)]
errorsMax_nb = [maximum([errorsMax_t_nb[idx_t][idx_f] for idx_t = eachindex(time_vec)]) for idx_f = 1:nPairs]


errorsL2 = [errorsL2_nb,errorsL2_b]
errorsMax = [errorsMax_nb, errorsMax_b]

fig2 = Figure(size=(700,350))
ax21 = Axis(fig2[1,1],xlabel="Δχ",ylabel="error",xscale=log10,yscale=log10,title="L2 Error")
for idx_b = 1:nBaths
    lines!(ax21,Dχs_b[pairsFirst],errorsL2[idx_b],label=baths[idx_b])
end
axislegend(ax21,position=:rb)
ax22 = Axis(fig2[1,2],xlabel="Δχ",ylabel="error",xscale=log10,yscale=log10,title="Maximum Error")
for idx_b = 1:nBaths
    lines!(ax22,Dχs_b[pairsFirst],errorsMax[idx_b],label=baths[idx_b])
end
axislegend(ax22,position=:rb)
fig2
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\etaConvErr.svg",fig2)


eocsL2_eta = [computeEOC(errorsL2[idx_b],Dχs[idx_b]) for idx_b = 1:nBaths]
eocsMax_eta = [computeEOC(errorsMax[idx_b],Dχs[idx_b]) for idx_b = 1:nBaths]

print("### ETA ###\nEOCs in L2 and max norm (idx_b x idx_f):\n")
@show eocsL2_eta
@show eocsMax_eta
avgL2 = sum.([eocsL2_eta[idx_b] for idx_b=1:nBaths])/(nFacs-2)
avgMax = sum.([eocsMax_eta[idx_b] for idx_b=1:nBaths])/(nFacs-2)
print("average:\nL2: $(avgL2),  Max: $(avgMax)\n")