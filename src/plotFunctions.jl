
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
function plotSensorData(numSensor::Sensors,time_vec::Vector{Float64})
    numTimeInds=nothing
    numLabel="eta_euler"
    title=""
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
            lines!(axs[i],time_vec,eta_sensor_real_interp[i].(time_vec),label="measurement")
        end
        if plotSWE
            lines!(axs[i],time_SWE,eta_sensor_SWE[i,:],label="eta_SWE")
        end
        lines!(axs[i],time_vec[numTimeInds],numSensor.data[i][numTimeInds],label=numLabel)
        xlims!(axs[i],(0,17))
        axislegend(axs[i],position=:lb)
    end
    Label(fig[begin-1,1:2], title, font = "Nimbus Sans Bold")

    return fig
end

function plotSensorData!(fig::GLMakie.Figure,numSensor::Sensors,time_vec::Vector{Float64};numTimeInds::Vector{Int}=nothing,numLabel::String="eta_euler")
    if numTimeInds === nothing
        numTimeInds = 1:length(time_vec)
    end
    #@show contents(fig) 
    for leg in contents(fig) 
        if  leg isa Legend
            delete!(leg)
        end
    end
    axs = []
    for content in contents(fig)
        if content isa Axis
            push!(axs,content)
        end
    end
    #@show axs
    for i =1:4
        lines!(axs[i],time_vec[numTimeInds],numSensor.data[i][numTimeInds],label=numLabel)
    end
    axislegend(axs[i],position=:lb)
end

function plotSurfaceOverTime(eta::Vector{Vector{Float64}},time_vec::Vector{Float64},χs::Vector{Float64},label::String)
    fig = Figure()
    ax = Axis(fig[1,1],xlabel="χ",ylabel="η",title="η on surface")
    t_sl = GLMakie.Slider(fig[2,1], range=eachindex(time_vec),startvalue=1)
    eta_t = GLMakie.lift(t_sl.value) do t
        eta[t]
    end
    M = maximum(maximum([abs.(eta[i]) for i = 1:length(eta)]))
    lines!(ax,χs,eta_t,label=label)
    ylims!(ax,-1.4*M,1.4*M)
    GLMakie.axislegend(ax,position=:lt)
    return fig
end

function plotSurfaceOverTime!(fig::GLMakie.Figure,eta::Vector{Vector{Float64}},time_vec::Vector{Float64},χs::Vector{Float64},label::String)
    for leg in fig.content 
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
    ax = axs[1]
    sls = []
    for content in fig.content
        if content isa GLMakie.Slider
            push!(sls,content)
        end
    end
    t_sl = sls[1]
    eta_t = GLMakie.lift(t_sl.value) do t
        eta[t]
    end
    lines!(ax,χs,eta_t,label=label)
    GLMakie.axislegend(ax,position=:lt);
end