using JLD2
filename = "NeumannUndamped.jld"
d = load(filename)
all_etas = d["all_etas"]
all_phis = d["all_phis"]