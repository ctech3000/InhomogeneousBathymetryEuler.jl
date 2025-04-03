using CairoMakie, JLD2, HDF5
using InhomogeneousBathymetryEuler
saveFigs = true

data = load("examples/results/plottingScripts/EulerDataSim_withBathyData.jld2")
sensorss = data["sensorss"]
time_vec = data["time_vec"]
domains = data["domains"]
tau05 = 17.45
t_ind_end_first_period = argmin(t_ind -> abs(tau05-time_vec[t_ind]),eachindex(time_vec))
first_period = time_vec[1:t_ind_end_first_period]
first_period_inds = collect(1:t_ind_end_first_period)
t_ind_start_second_period = argmin(t_ind -> abs(2*tau05-time_vec[t_ind]),eachindex(time_vec))
t_ind_end_second_period = argmin(t_ind -> abs(3*tau05-time_vec[t_ind]),eachindex(time_vec))
second_period = time_vec[t_ind_start_second_period:t_ind_end_second_period]
second_period_inds = collect(t_ind_start_second_period:t_ind_end_second_period)
fiveSecInd = 401

data_nb = load("examples/results/plottingScripts/EulerDataSim_noBathyData.jld2")
sensorss_nb = data_nb["sensorss"]

fig1_ = plotSensorData(sensorss_nb[3],first_period,numLabel="no bath.",bathy=false,plotSWE=false,plotExp=false,numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
plotSensorData!(fig1_,sensorss[1,3],first_period,numLabel="with bath.",numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
Legend(fig1_[3,1:2],fig1_.content[1],orientation=:horizontal)
resize!(fig1_.scene,(750,500))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\modelComp_withAndWithoutBathy.svg",fig1_)
end
fig1_

fig1 = plotSensorData(sensorss[1,2],first_period,numLabel=L"Euler, $w_{nb}$",bathy=true,plotSWE=false,numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
plotSensorData!(fig1,sensorss[2,2],first_period,numLabel=L"Euler, $w_b$",numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
Legend(fig1[3,1:2],fig1.content[1],orientation=:horizontal)
resize!(fig1.scene,(750,550))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\modelComp_withBathy_sensorWave.svg",fig1)
end
fig1

fig2 = plotSensorData(sensorss[1,1],first_period,numLabel=L"direct$$",bathy=true,plotSWE=false,numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
plotSensorData!(fig2,sensorss[1,2],first_period,numLabel=L"relaxed, $L_D=2λ$",numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
plotSensorData!(fig2,sensorss[1,3],first_period,numLabel=L"relaxed, $L_D=4λ$",numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
Legend(fig2[3,1:2],fig2.content[1],orientation=:horizontal)
resize!(fig2.scene,(750,550))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\modelComp_withBathy_waveGeneration.svg",fig2)
end
fig2

fig3 = plotSensorData(sensorss[1,3],first_period,numLabel="Euler",bathy=true,plotSWE=true,numTimeInds=collect(fiveSecInd:t_ind_end_first_period))
Legend(fig3[3,1:2],fig3.content[1],orientation=:horizontal)
resize!(fig3.scene,(750,550))
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\modelComp_withBathy_withSWE.svg",fig3)
end
fig3

