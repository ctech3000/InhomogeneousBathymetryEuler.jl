using CairoMakie, JLD2
using InhomogeneousBathymetryEuler

data = load("examples/results/plottingScripts/convManSolData.jld2")
errorsL2 = data["errorsL2"]
errorsMax = data["errorsMax"]
Dχs = data["Dχs"]
nGrids, nBaths = size(errorsL2)
baths = ["b₀","b₁","b₂","b₃"]
nBaths = size(errorsL2,2)
nFacs = size(errorsL2,1)


eocsL2 = [computeEOC(errorsL2[:,idx_b],Dχs[:,idx_b]) for idx_b = 1:nBaths]
eocsMax = [computeEOC(errorsMax[:,idx_b],Dχs[:,idx_b]) for idx_b = 1:nBaths]

CairoMakie.activate!()
set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))

fig1 = Figure(size=(700,300))
ax11 = Axis(fig1[1,1],xlabel="Δχ",ylabel="error",xscale=log10,yscale=log10,title="L2 Error")
for idx_b = 1:nBaths
    lines!(ax11,Dχs[:,idx_b],errorsL2[:,idx_b],label=baths[idx_b])
end
axislegend(ax11,position=:rb)
ax12 = Axis(fig1[1,2],xlabel="Δχ",ylabel="error",xscale=log10,yscale=log10,title="Maximum Error")
for idx_b = 1:nBaths
    lines!(ax12,Dχs[:,idx_b],errorsMax[:,idx_b],label=baths[idx_b])
end
axislegend(ax12,position=:rb)

avgL2 = sum.([eocsL2[idx_b] for idx_b=1:nBaths])/(nFacs-1)
avgMax = sum.([eocsMax[idx_b] for idx_b=1:nBaths])/(nFacs-1)

save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\manSolErr.svg",fig1)
fig1