using CairoMakie, JLD2, HDF5
using InhomogeneousBathymetryEuler
CairoMakie.activate!()
saveFigs = false

data = load("examples/results/presentationData/dataComp_flat.jld2")
sensors_flat = data["sensors"]
time_vec = data["time_vec"]
eta_flat = data["all_etas"]
tau05 = 17.45
t_ind_end_first_period = argmin(t_ind -> abs(tau05-time_vec[t_ind]),eachindex(time_vec))
first_period = time_vec[1:t_ind_end_first_period]

fiveSecInd = 201

fig1 = plotSensorData_alt(sensors_flat,first_period,numLabel="Euler",bathy=false,numTimeInds=collect(fiveSecInd:t_ind_end_first_period),sweLinestyle=:solid)
Legend(fig1[3,1:2],fig1.content[1],orientation=:horizontal)
resize!(fig1.scene,(1050,550))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-presentation\\images\\compFlat.svg",fig1)
end
xlims!(fig1.content[3],(13,15))
ylims!(fig1.content[3],(0.0025,0.007))
xlims!(fig1.content[4],(12.5,15))
ylims!(fig1.content[4],(-0.006,-0.001))
fig1
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-presentation\\images\\compFlat_close.svg",fig1)
end


data = load("examples/results/presentationData/dataComp_gauss.jld2")
sensors_gauss = data["sensors"]
eta_gauss = data["all_etas"]

tau05 = 15
t_ind_end_first_period = argmin(t_ind -> abs(tau05-time_vec[t_ind]),eachindex(time_vec))
first_period = time_vec[1:t_ind_end_first_period]

saveFigs = true
fig2 = plotSensorData_alt(sensors_gauss,first_period,numLabel="Euler",bathy=true,numTimeInds=collect(fiveSecInd:t_ind_end_first_period),sweLinestyle=:solid)
Legend(fig2[3,1:2],fig2.content[1],orientation=:horizontal)
resize!(fig2.scene,(1050,550))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-presentation\\images\\compGauss.svg",fig2)
end
fig2

xlims!(fig2.content[3],(11,12))
ylims!(fig2.content[3],(0.0025,0.0051))
xlims!(fig2.content[4],(12,13.2))
ylims!(fig2.content[4],(0.0025,0.0051))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-presentation\\images\\compGauss_close.svg",fig2)
end