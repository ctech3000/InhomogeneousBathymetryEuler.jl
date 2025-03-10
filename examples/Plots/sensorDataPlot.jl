using CairoMakie, DelimitedFiles

CairoMakie.activate!(type="svg")

filename = "C:\\Users\\chris\\.julia\\dev\\InhomogeneousBathymetryEuler.jl\\examples\\Data_Sensors\\Without_Bathymetry\\Heat2.txt"
eta_sensor = Vector{Vector{Float64}}(undef,4)
for sensor_idx = 1:4
    eta_sensor[sensor_idx] = convert(Vector{Float64},readdlm(filename)[2:end,sensor_idx])/100
end
time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])

set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))

time_interval = 2800:5:5200
fig = Figure(size=(400,230))
ax = Axis(fig[1,1],xlabel="time (s)",ylabel="displacement (m)")
for sensor = [1,4]
    lines!(ax,time_sensor[time_interval],eta_sensor[sensor][time_interval],label="Sensor $(sensor)")
end
axislegend(position=:lt)
fig
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\sensorData.svg",fig)