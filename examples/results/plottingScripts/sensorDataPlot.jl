# produces figures 2
using CairoMakie, DelimitedFiles

CairoMakie.activate!(type="svg")
virColors=[1,9]
MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))

eta_sensor_real = [zeros(Float64,10000) for sensor = 1:4]
time_sensor = Vector{Float64}(undef,10000)
bathy = false
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

time_interval = 2800:5:5200
fig = Figure(size=(400,230),font="CMU")
ax = Axis(fig[1,1],xlabel="time /s",ylabel="displacement /m")
for (i,sensor) = enumerate([1,4])
    lines!(ax,time_sensor[time_interval],eta_sensor_real[sensor][time_interval],label="Sensor $(sensor)",color=virColors[i], colormap=:viridis, colorrange=(1,10))
end
axislegend(position=:lt)
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\sensorData.svg",fig)
fig