using InhomogeneousBathymetryEuler, JLD2

irreg = true

if irreg
    filename = "examples/Plots/dampingLayerDataIrregEuler.jld2"
else
    filename = "examples/Plots/dampingLayerDataRegEuler.jld2"
end
d = load(filename)

