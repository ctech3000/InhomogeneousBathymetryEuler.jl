using JLD
filename = "DirichletDamped.jld"
d = load(filename)
all_etas = d["all_etas"]
all_phis = d["all_phis"]