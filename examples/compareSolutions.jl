using DelimitedFiles, Interpolations

function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I
end


filename = "C:\\Users\\chris\\.julia\\dev\\InhomogeneousBathymetryEuler.jl\\examples\\Data_Sensors\\Without_Bathymetry\\Heat1.txt"
eta_sensor = [convert(Vector{Float64},readdlm(filename)[2:end,sensor])/100 for sensor = 1:4]
time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])
eta_sensor_interp = [linear_interpolation(time_sensor,eta_sensor[sensor], extrapolation_bc = Flat()) for sensor = 1:4]

eta_wave = [[-1/g*analyticPotential_dt(sensors.sensors_pos_x[sensor],0.0,t,0.3,wave) for t in time_vec] for sensor = 1:4]

eta_num = sensors.data

fig1 = Figure()
axs1 = [Axis(fig1[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs1[i],time_vec,eta_wave[i],label="eta_wave")
    lines!(axs1[i],time_vec,eta_num[i],label="eta_num")
    axislegend(axs1[i],position=:lb)
end
Label(fig1[begin-1,1:2], "Comparison Computed/Analytic, reg.", font = "Nimbus Sans Bold")
fig1

fig2 = Figure()
axs2 = [Axis(fig2[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
t_sl = Slider(fig2[3, 1:2], range = 0:0.2:30, startvalue = 0)
for i =1:4
    eta_sensor_sl = lift(t_sl.value) do t
        eta_sensor_interp[i].(time_vec.+t)
    end
    lines!(axs2[i],time_vec,eta_sensor_sl,label="sensor_eta")
    lines!(axs2[i],time_vec,eta_num[i],label="eta_num")
    axislegend(axs2[i],position=:lb)
end
Label(fig2[begin-1,1:2], "Comparison Computed/Sensor, reg.", font = "Nimbus Sans Bold")
fig2