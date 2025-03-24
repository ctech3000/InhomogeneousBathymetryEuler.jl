using CairoMakie, JLD2
using InhomogeneousBathymetryEuler

data = load("examples/results/plottingScripts/convManSolData.jld2")
errorsL2 = data["errorsL2"]
errorsMax = data["errorsMax"]
Dχs = data["Dχs"]
nGrids, nBaths = size(errorsL2)
baths = ["Flat","s.G.","b.G.","s.S.","b.S."]


eocsL2 = [computeEOC(errorsL2[:,idx_b],Dχs[:,idx_b]) for idx_b = 1:nBaths]
eocsMax = [computeEOC(errorsMax[:,idx_b],Dχs[:,idx_b]) for idx_b = 1:nBaths]

set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))

fig1 = Figure(size=(700,400))
ax11 = Axis(fig1[1,1],xlabel="Δχ",ylabel="error",xscale=log2,yscale=log2,title="L2 error")
for idx_b = 1:nBaths
    lines!(ax11,Dχs[:,idx_b],errorsL2[:,idx_b],label=baths[idx_b])
end
axislegend(ax11,position=:rb)
ax12 = Axis(fig1[1,2],xlabel="Δχ",ylabel="error",xscale=log2,yscale=log2,title="max error")
for idx_b = 1:nBaths
    lines!(ax12,Dχs[:,idx_b],errorsMax[:,idx_b],label=baths[idx_b])
end
axislegend(ax12,position=:rb)
fig1