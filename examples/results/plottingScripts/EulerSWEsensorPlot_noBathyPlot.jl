using CairoMakie, JLD2, HDF5
using InhomogeneousBathymetryEuler
saveFigs = false

data = load("examples/results/plottingScripts/EulerDataSim_noBathyData.jld2")
sensorss = data["sensorss"]
time_vec = data["time_vec"]
tau05 = 17.45
t_ind_end_first_period = argmin(t_ind -> abs(tau05-time_vec[t_ind]),eachindex(time_vec))
first_period = time_vec[1:t_ind_end_first_period]
first_period_inds = collect(1:t_ind_end_first_period)
t_ind_start_second_period = argmin(t_ind -> abs(2*tau05-time_vec[t_ind]),eachindex(time_vec))
t_ind_end_second_period = argmin(t_ind -> abs(3*tau05-time_vec[t_ind]),eachindex(time_vec))
second_period = time_vec[t_ind_start_second_period:t_ind_end_second_period]
second_period_inds = collect(t_ind_start_second_period:t_ind_end_second_period)
fiveSecInd = 401

fig1 = plotSensorData(sensorss[2],first_period,numLabel="Euler",bathy=false,numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
Legend(fig1[3,1:2],fig1.content[1],orientation=:horizontal)
resize!(fig1.scene,(750,550))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\modelComp_noBathy.svg",fig1)
end
fig1