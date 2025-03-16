using GLMakie, DelimitedFiles, Interpolations
GLMakie.activate!()
time_ind = 1:length(time_vec)


###
#= plot sensor data with bathy (true) or without (false) =#
bathy = false
#= Data_Sensors folder with sensor data has to be in examples folder  =#
###

# plot surface displacement
amp = 0.005
f1 = Figure()
ax1 = Axis(f1[1,1],xlabel="χ",ylabel="η",title = "η on whole surface")

t_sl = Slider(f1[2, 1], range = time_ind, startvalue = 1)
eta = lift(t_sl.value) do t
    all_etas[t]
end
lines!(χs,eta)
ylims!(ax1,-1.5*amp,1.5*amp)
f1

# plot phi on surface
f2 = Figure()
g = 9.81
phi_amp = g*amp/getFreq(wave)
ax2 = Axis(f2[1,1],xlabel="χ",ylabel="ϕ",title="ϕ on whole surface")
time_ind = 1:length(time_vec)
t_sl = Slider(f2[2, 1], range = time_ind, startvalue = 1)
phi = lift(t_sl.value) do t
    all_phis[t][:,1]
end
lines!(χs,phi)
ylims!(ax2,-1.5*phi_amp,1.5*phi_amp)
f2

# plot simulated and real sensor data
function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I
end
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
eta_sensor_real_interp = [linear_interpolation(time_sensor.-27,eta_sensor_real[sensor], extrapolation_bc = Flat()) for sensor = 1:4]
f3 = Figure(size=(1280,720))
axs3 = [Axis(f3[linToCart(i,2,2)...],xlabel="t", ylabel="η",title="Sensor $(i)") for i = 1:4]
for i =1:4
    lines!(axs3[i],time_vec,eta_sensor_real_interp[i].(time_vec),label="measurement")
    lines!(axs3[i],time_vec,sensors.data[i],label="eta_num")
    axislegend(axs3[i],position=:lb)
end
if bathy
    out = ""
else
    out = "out"
end
Label(f3[begin-1,1:2], "Comparison num. Euler/Sensor, with"*out*" bathymetry", font = "Nimbus Sans Bold")
f3
