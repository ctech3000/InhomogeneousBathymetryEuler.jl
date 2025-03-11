using InhomogeneousBathymetryEuler, JLD2, HDF5, DelimitedFiles, Interpolations, CairoMakie

function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I
end

done_once = false
if !done_once
    #read swe data
    M = 100
    fid = h5open("examples/Plots/nobathyHeat_mean_kappa2e-01_T=20_M=100.hdf5", "r")
    eta_sensor234_SWE = read(fid,"H_sensor")[:,1:15001].-0.3
    eta_sensor1_SWE = read(fid,"h")[1,1:15001].-0.3
    eta_sensor_SWE = vcat(eta_sensor1_SWE',eta_sensor234_SWE)
    time_SWE = read(fid,"t_array")[1:15001].+3.02

    #read sensor data
    eta_sensor_real = [zeros(Float64,10000) for sensor = 1:4]
    for heat_idx = 1:20
        filename = "examples/Data_Sensors/Without_Bathymetry/Heat$(heat_idx).txt"
        for sensor_idx = 1:4
            global eta_sensor_real[sensor_idx] += convert(Vector{Float64},readdlm(filename)[2:end,sensor_idx])/100
        end
        global time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])
    end
    for sensor_idx = 1:4
        eta_sensor_real[sensor_idx] = eta_sensor_real[sensor_idx]/20
    end
    eta_sensor_real_interp = [linear_interpolation(time_sensor.-27,eta_sensor_real[sensor], extrapolation_bc = Flat()) for sensor = 1:4]

    #read euler data
    euler_time_ind = 1:301
    d = JLD2.load("examples/Plots/compSolData.jld2")
    sensors_reg = d["sensors_reg"]
    sensors_irreg = d["sensors_irreg"]
    time_vec = d["time_vec"]
    wave_reg = d["wave_reg"]
    wave_irreg = d["wave_irreg"]
    eta_sensor_euler_reg = [sensors_reg.data[sensor_idx][euler_time_ind] for sensor_idx = 1:4]
    eta_sensor_euler_irreg = [sensors_irreg.data[sensor_idx][euler_time_ind] for sensor_idx = 1:4]
    time_vec = time_vec[euler_time_ind]

    #compute analytic TrueSolution
    eta_wave_reg = [[-1/GRAV*analyticPotential_dt(sensors_reg.sensors_pos_x[sensor],0.0,t,0.3,wave_reg) for t in time_vec] for sensor = 1:4]
    eta_wave_irreg = [[-1/GRAV*analyticPotential_dt(sensors_irreg.sensors_pos_x[sensor],0.0,t,0.3,wave_irreg) for t in time_vec] for sensor = 1:4]
    done_once = true
end

fig1 = Figure()
axs1 = [Axis(fig1[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs1[i],time_vec,eta_wave_reg[i],label="eta_wave_reg")
    lines!(axs1[i],time_vec,eta_sensor_euler_reg[i],label="eta_num_reg")
    axislegend(axs1[i],position=:lb)
end
Label(fig1[begin-1,1:2], "Comparison num. Euler/Ana. Euler, reg.", font = "Nimbus Sans Bold")
fig1

fig2 = Figure()
axs2 = [Axis(fig2[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs2[i],time_vec,eta_wave_irreg[i],label="eta_wave_irreg")
    lines!(axs2[i],time_vec,eta_sensor_euler_irreg[i],label="eta_num_irreg")
    axislegend(axs2[i],position=:lb)
end
Label(fig2[begin-1,1:2], "Comparison num. Euler/Ana. Euler, irreg.", font = "Nimbus Sans Bold")
fig2

fig3 = Figure()
axs3 = [Axis(fig3[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs3[i],time_vec,eta_sensor_real_interp[i].(time_vec),label="measurement")
    lines!(axs3[i],time_SWE,eta_sensor_SWE[i,:],label="eta_SWE")
    lines!(axs3[i],time_vec.+3.2,eta_sensor_euler_reg[i],label="eta_num_reg")
    axislegend(axs3[i],position=:lb)
end
Label(fig3[begin-1,1:2], "Comparison num. Euler, reg./SWE", font = "Nimbus Sans Bold")
fig3

fig4 = Figure(size=(1280,720))
axs4 = [Axis(fig4[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs4[i],time_vec,eta_sensor_real_interp[i].(time_vec),label="measurement")
    lines!(axs4[i],time_SWE,eta_sensor_SWE[i,:],label="eta_SWE")
    lines!(axs4[i],time_vec,eta_sensor_euler_irreg[i],label="eta_num_irreg")
    axislegend(axs4[i],position=:lb)
end
Label(fig4[begin-1,1:2], "Comparison num. Euler,irreg./SWE, without bath", font = "Nimbus Sans Bold")
fig4