using InhomogeneousBathymetryEuler
using JLD2, GLMakie, HDF5, DelimitedFiles, Interpolations

function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I
end

function loadSWE(bathy::Bool)
    if bathy
        filename="examples/Plots/sim_data_Heat_meanBathy_ExactRamp_T=20_M=100.hdf5"
    else
        filename="examples/Plots/nobathyHeat_mean_kappa2e-01_T=20_M=100.hdf5"
    end
    fid = h5open(filename, "r")
    eta_sensor234_SWE = read(fid,"H_sensor")[:,1:15001].-0.3
    eta_sensor1_SWE = read(fid,"h")[1,1:15001].-0.3
    eta_sensor_SWE = vcat(eta_sensor1_SWE',eta_sensor234_SWE)
    time_SWE = read(fid,"t_array")[1:15001].+3.02
    return eta_sensor_SWE, time_SWE
end

function loadRealSensorData(bathy::Bool)
    eta_sensor_real = [zeros(Float64,10000) for sensor = 1:4]
    for heat_idx = 1:20
        if bathy
            filename = "examples/Data_Sensors/With_Bathymetry/Heat$(heat_idx).txt"
        else
            filename = "examples/Data_Sensors/Without_Bathymetry/Heat$(heat_idx).txt"
        end
        for sensor_idx = 1:4
            global eta_sensor_real[sensor_idx] += convert(Vector{Float64},readdlm(filename)[2:end,sensor_idx])/100
        end
        global time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])
    end
    for sensor_idx = 1:4
        eta_sensor_real[sensor_idx] = eta_sensor_real[sensor_idx]/20
    end
    eta_sensor_real_interp = [Interpolations.linear_interpolation(time_sensor.-27,eta_sensor_real[sensor], extrapolation_bc = Interpolations.Flat()) for sensor = 1:4]
    return eta_sensor_real_interp
end
#numTimeInds::Vector{Int}=nothing,numLabel::String="eta_euler",title::String="",plotSWE::Bool=true,plotRealSensor::Bool=true,bathy::Bool=true
function plotSensorData(numSensor::Sensors,time_vec::Vector{Float64},numLabel::String,title::String)
    numTimeInds=nothing
    #numLabel="eta_euler"
    #title=""
    plotSWE=true
    plotRealSensor=true
    bathy=true
    if plotSWE
        eta_sensor_SWE, time_SWE = loadSWE(bathy)
    end
    if plotRealSensor
        eta_sensor_real_interp = loadRealSensorData(bathy)
    end
    if numTimeInds === nothing
        numTimeInds = 1:length(time_vec)
    end
    fig = Figure()
    axs = [Axis(fig[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
    for i =1:4
        if plotRealSensor
            lines!(axs[i],time_vec,eta_sensor_real_interp[i].(time_vec),label="measurement",linewidth=2.5)
        end
        if plotSWE
            lines!(axs[i],time_SWE,eta_sensor_SWE[i,:],label="eta_SWE",linestyle=:dash,color=:gray)
        end
        lines!(axs[i],time_vec[numTimeInds],numSensor.data[i][numTimeInds],label=numLabel)
        xlims!(axs[i],(0,17))
        axislegend(axs[i],position=:lb)
    end
    Label(fig[begin-1,1:2], title, font = "Nimbus Sans Bold")

    return fig
end
#;numTimeInds::Vector{Int}=nothing,numLabel::String="eta_euler"
function plotSensorData!(fig::GLMakie.Figure,numSensor::Sensors,time_vec::Vector{Float64},numLabel::String)
    numTimeInds=nothing
    #numLabel="eta_euler"
    if numTimeInds === nothing
        numTimeInds = 1:length(time_vec)
    end
    for leg in reverse(fig.content)
        if  leg isa Legend
            delete!(leg)
        end
    end
    axs = []
    for content in fig.content
        if content isa Axis
            push!(axs,content)
        end
    end
    for i =1:4
        lines!(axs[i],time_vec[numTimeInds],numSensor.data[i][numTimeInds],label=numLabel)
        xlims!(axs[i],(0,17))
        axislegend(axs[i],position=:lb)
    end
end

d = JLD2.load("examples/Plots/relaxationLengthData.jld2")
num_sensors_heats = d["sensors_heats"]
time_vec = d["time_vec"]
extension_vec = d["extension_vec"]
wave,_,_ = IrregWave("examples/irregWaveData_noBathy.jld2",inflowDepth=0.3) # irreg wave
nHeats = length(num_sensors_heats)
eta_num_sensor_heats = [[num_sensors_heats[heat].data[sensor_idx] for sensor_idx = 1:4] for heat = 1:nHeats]

fig = plotSensorData(num_sensors_heats[1],time_vec,"euler, L_RX=$(extension_vec[1])λ","400 points per λ")
for i = [3,5,7]
    plotSensorData!(fig,num_sensors_heats[i],time_vec,"euler, L_RX=$(extension_vec[i])λ")
end
fig

d = JLD2.load("examples/Plots/relaxationLengthDataFiner.jld2")
num_sensors_heats = d["sensors_heats"]
time_vec = d["time_vec"]
extension_vec = d["extension_vec"]
wave,_,_ = IrregWave("examples/irregWaveData_noBathy.jld2",inflowDepth=0.3) # irreg wave
nHeats = length(num_sensors_heats)
eta_num_sensor_heats = [[num_sensors_heats[heat].data[sensor_idx] for sensor_idx = 1:4] for heat = 1:nHeats]

fig2 = plotSensorData(num_sensors_heats[1],time_vec,"euler, L_RX=$(extension_vec[1])λ","800 points per λ")
for i = [3,5,7]
    plotSensorData!(fig2,num_sensors_heats[i],time_vec,"euler, L_RX=$(extension_vec[i])λ")
end
fig2