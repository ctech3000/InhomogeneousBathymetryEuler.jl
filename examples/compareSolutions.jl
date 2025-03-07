using DelimitedFiles, Interpolations, GLMakie

function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I
end

eta_sensor = [zeros(Float64,10000) for sensor = 1:4]
for heat_idx = 1:20
    filename = "C:\\Users\\chris\\.julia\\dev\\InhomogeneousBathymetryEuler.jl\\examples\\Data_Sensors\\Without_Bathymetry\\Heat$(heat_idx).txt"
    for sensor_idx = 1:4
        global eta_sensor[sensor_idx] += convert(Vector{Float64},readdlm(filename)[2:end,sensor_idx])/100
    end
    global time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])
end
for sensor_idx = 1:4
    eta_sensor[sensor_idx] = eta_sensor[sensor_idx]/20 .- amp0
end
eta_sensor_interp = [linear_interpolation(time_sensor,eta_sensor[sensor], extrapolation_bc = Flat()) for sensor = 1:4]

eta_wave_reg = [[-1/GRAV*analyticPotential_dt(sensors.sensors_pos_x[sensor],0.0,t,0.3,wave_reg) for t in time_vec] for sensor = 1:4]
eta_wave_irreg = [[-1/GRAV*analyticPotential_dt(sensors.sensors_pos_x[sensor],0.0,t,0.3,wave_irreg) for t in time_vec] for sensor = 1:4]

eta_num_reg = sensors_reg.data
eta_num_irreg = sensors_irreg.data

fig1 = Figure()
axs1 = [Axis(fig1[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs1[i],time_vec,eta_wave_reg[i],label="eta_wave_reg")
    lines!(axs1[i],time_vec,eta_num_reg[i],label="eta_num_reg")
    axislegend(axs1[i],position=:lb)
end
Label(fig1[begin-1,1:2], "Comparison Computed/Analytic, reg.", font = "Nimbus Sans Bold")
fig1

fig2 = Figure()
axs2 = [Axis(fig2[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
t_sl = Slider(fig2[3, 1:2], range = 0:0.05:50, startvalue = 0)
for i =1:4
    eta_sensor_sl = lift(t_sl.value) do t
        eta_sensor_interp[i].(time_vec.+t)
    end
    lines!(axs2[i],time_vec,eta_sensor_sl,label="sensor_eta")
    lines!(axs2[i],time_vec,eta_num_reg[i],label="eta_num_reg")
    axislegend(axs2[i],position=:lb)
end
Label(fig2[begin-1,1:2], "Comparison Computed/Sensor, reg.", font = "Nimbus Sans Bold")
fig2

fig3 = Figure()
axs3 = [Axis(fig3[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs3[i],time_vec,eta_wave_irreg[i],label="eta_wave_irreg")
    lines!(axs3[i],time_vec,eta_num_irreg[i],label="eta_num_irreg")
    axislegend(axs3[i],position=:lb)
end
Label(fig3[begin-1,1:2], "Comparison Computed/Analytic, irreg.", font = "Nimbus Sans Bold")
fig3

fig3_5 = Figure()
axs3_5 = [Axis(fig3_5[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
t_sl = Slider(fig3_5[3, 1:2], range = 0:0.05:50, startvalue = 0)
for i =1:4
    eta_sensor_sl = lift(t_sl.value) do t
        eta_sensor_interp[i].(time_vec.+t)
    end
    lines!(axs3_5[i],time_vec,eta_sensor_sl,label="sensor_eta")
    lines!(axs3_5[i],time_vec,eta_wave_irreg[i],label="eta_wave_irreg")
    axislegend(axs3_5[i],position=:lb)
end
Label(fig3_5[begin-1,1:2], "Comparison Sensor/Analytic, irreg.", font = "Nimbus Sans Bold")
fig3_5

fig4 = Figure()
axs4 = [Axis(fig4[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
t_sl = Slider(fig4[3, 1:2], range = 0:0.05:50, startvalue = 0)
for i =1:4
    eta_sensor_sl = lift(t_sl.value) do t
        eta_sensor_interp[i].(time_vec.+t)
    end
    lines!(axs4[i],time_vec,eta_sensor_sl,label="sensor_eta")
    lines!(axs4[i],time_vec,eta_num_irreg[i],label="eta_num_irreg")
    axislegend(axs4[i],position=:lb)
end
Label(fig4[begin-1,1:2], "Comparison Computed/Sensor, irreg.", font = "Nimbus Sans Bold")
fig4

fig5 = Figure()
axs5 = [Axis(fig5[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
t_sl = Slider(fig5[3, 1:2], range = 0:0.05:50, startvalue = 0)
for i =1:4
    eta_sensor_sl = lift(t_sl.value) do t
        eta_sensor_interp[i].(time_vec.+t)
    end
    lines!(axs5[i],time_vec,eta_sensor_sl,label="sensor_eta")
    lines!(axs5[i],time_vec,eta_num_reg[i],label="eta_num_reg")
    lines!(axs5[i],time_vec,eta_num_irreg[i],label="eta_num_irreg")
    axislegend(axs5[i],position=:lb)
end
Label(fig5[begin-1,1:2], "Comparison Computed/Sensor, irreg. and reg.", font = "Nimbus Sans Bold")
fig5