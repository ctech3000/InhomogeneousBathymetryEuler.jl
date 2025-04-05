using CairoMakie, DelimitedFiles

CairoMakie.activate!(type="svg")
virColors=[1,9]
MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))

filename = "C:\\Users\\chris\\.julia\\dev\\InhomogeneousBathymetryEuler.jl\\examples\\Data_Sensors\\Without_Bathymetry\\Heat2.txt"
eta_sensor = Vector{Vector{Float64}}(undef,4)
for sensor_idx = 1:4
    eta_sensor[sensor_idx] = convert(Vector{Float64},readdlm(filename)[2:end,sensor_idx])/100
end
time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])

time_interval = 2800:5:5200
fig = Figure(size=(400,230),font="CMU")
ax = Axis(fig[1,1],xlabel="time /s",ylabel="displacement /m")
for (i,sensor) = enumerate([1,4])
    lines!(ax,time_sensor[time_interval],eta_sensor[sensor][time_interval],label="Sensor $(sensor)",color=virColors[i], colormap=:viridis, colorrange=(1,10))
end
axislegend(position=:lt)
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\sensorData.svg",fig)
fig