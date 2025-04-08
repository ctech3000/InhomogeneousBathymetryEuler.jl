### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ac91d146-a4d4-44a7-ace7-fdb770f0b367
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate(Base.current_project())
    using Revise # Revise must be in your LOAD_PATH
end

# ╔═╡ 7e8dab8d-fc6e-4a67-8ba1-87e48955b627
using InhomogeneousBathymetryEuler, PlutoUI, GLMakie, FFTW, DelimitedFiles, JLD2

# ╔═╡ 2deecc02-48f9-46d3-8ac2-fba747562cbb
html"""
<style>
input[type*="range"] {
	width: 100%;
}
</style>
"""

# ╔═╡ 1834c732-00e8-4cbf-9869-7853f9174af3
md"""
Select with or without bathymetry.
"""

# ╔═╡ 58d9e3eb-2187-4eee-8ea1-457429b38fda
@bind bath Select(["with bathymetry", "without bathymetry"])

# ╔═╡ cad908ae-1fa3-4a04-ae91-282d34ef4562
md"""
Select number of heats to be averaged.
"""

# ╔═╡ 902d2a0a-cd9c-4f56-a04e-3c155fcea3e4
@bind nHeats PlutoUI.Slider(1:20,default = 20,show_value=true)

# ╔═╡ 26daba59-0f1d-4dbe-8e70-1cc78ddb2f5e
begin
	eta = zeros(Float64,10000)
	time_sensor = zeros(Float64,10000)
	for heat_idx = 1:nHeats
		if bath == "with bathymetry"
			filename = "Data_Sensors/With_Bathymetry/Heat$(heat_idx).txt"
		else
			filename = "Data_Sensors/Without_Bathymetry/Heat$(heat_idx).txt"
		end
		sensor = 1
		sensor_pos = [0.0,2.0,4.0,6.0]
		x = sensor_pos[sensor]
		z = 0.0
		global eta += convert(Vector{Float64},readdlm(filename)[2:end,sensor])/100
		global time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])
	end
	eta /= nHeats
end;

# ╔═╡ e4cc7114-3b89-48e7-847a-06bd961e10d8
md"""
Select start time.
"""

# ╔═╡ 66c4193d-fbc2-4a67-84dd-ef18652171b6
begin
	@bind t_start PlutoUI.Slider(time_sensor, default = 27.0, show_value=true)
end

# ╔═╡ f8c7bc9c-76ab-4c4c-aff6-238e622a21bf
md"""
Select end time.
"""

# ╔═╡ b5634eba-2bc9-4de5-a13e-afb4a7cb38e7
@bind t_end PlutoUI.Slider(time_sensor, default = 44.45, show_value=true)

# ╔═╡ 8ae9e601-bbdd-4e72-b8b0-13218d77e91c
begin
	t_start_ind = findall(x->x==t_start,time_sensor)[1]
	t_end_ind = findall(x->x==t_end,time_sensor)[1]
	selected_time_inds = t_start_ind:t_end_ind
	selected_time = time_sensor[selected_time_inds]
	selected_eta = eta[selected_time_inds]
end;

# ╔═╡ a0af295d-2850-41f6-a6bc-e5576b957bd7
md"""
Sensor input into Simulation:
"""

# ╔═╡ 5f652052-e23b-409d-8776-074e13a49699
begin
	N_t = length(selected_time)-1
	dt = (t_end - t_start)/(N_t+1)
	shifted_time = selected_time[1:end-1] .- selected_time[1] .+ 1/2*dt
	t_N = shifted_time[end]
	shifted_eta = [(selected_eta[i]+selected_eta[i+1])/2 for i = 1:N_t]
	coeff = dct(shifted_eta)
	omegas = [pi*k/t_N for k=0:N_t-1]
	amps = sqrt(2/(N_t))*coeff
	cos_reco = [sum(amps.*cos.(omegas*t)) for t in shifted_time]
end;

# ╔═╡ 0860d49f-33b1-4385-8aad-5d06c44715d4
md"""
Select maximum k (number of wave components).
"""

# ╔═╡ 3c85ca95-00a7-4307-ae5f-9c6f9833d2f9
@bind max_k PlutoUI.Slider(eachindex(coeff), default = 60, show_value=true)

# ╔═╡ d4cb9ea5-7b1e-47d6-ba81-d434c4624727
begin
	filtered_amps = amps[1:max_k]
	filtered_omegas = omegas[1:max_k]
	filtered_cos_reco = [sum(filtered_amps.*cos.(filtered_omegas*t)) for t in shifted_time]
end;

# ╔═╡ 4f9c2646-6552-4cfa-9ff8-b8eefba7de86
begin
	phases = zeros(size(filtered_amps))
	wave = IrregWave(filtered_amps[2:end],filtered_omegas[2:end],phases[2:end],hasFadeIn=false,inflowDepth=0.3)
	comp_eta = -1/9.81*[analyticPotential_dt(0,0,t,0.3,wave) for t in shifted_time]
	fig4 = Figure()
	ax4 = Axis(fig4[1,1])
	lines!(ax4,shifted_time,selected_eta[1:end-1],label="orig.")
	lines!(ax4,shifted_time,filtered_cos_reco, label="filtered cos_reco")
	lines!(ax4,shifted_time,comp_eta,label="comp_eta")
	axislegend(position = :lt)
	fig4
end

# ╔═╡ 4bc4bd94-17c3-4181-a70e-f02f9cad4f02
begin 
	using CairoMakie
	CairoMakie.activate!()
	set_theme!(fonts = (; regular = "Liberation Serif", bold = "Liberation Serif Bold"))
	irreg_wave = IrregWave(filtered_amps[2:end],filtered_omegas[2:end],phases[2:end],hasFadeIn=true,inflowDepth=0.3)
	irreg_eta = -1/9.81*[analyticPotential_dt(0,0,t,0.3,irreg_wave) for t in shifted_time]
end;

# ╔═╡ 1ecbc240-6615-4249-8032-f569f1fb8524
begin 
	lines(selected_time,selected_eta)
end

# ╔═╡ 3edc7820-70e5-4128-a472-a7ab7a1c095f
begin
	fig3 = Figure(size = (600,900))
	ax3 = Axis(fig3[1,1],xlabel="t",title="reconstruction")
	lines!(ax3,shifted_time,selected_eta[1:end-1],label="orig.")
	lines!(ax3,shifted_time,filtered_cos_reco, label="filtered cos_reco")
	axislegend(position = :lb)
	ax3_5 = Axis(fig3[2,1],xlabel="ω",ylabel="amp",title="coefficients")
	lines!(ax3_5,filtered_omegas,filtered_amps)
	fig3
end

# ╔═╡ 3ac0bd7a-0248-4e18-b64e-32ec0d9c7404
save_data = false

# ╔═╡ 7529d603-1211-464c-903c-e45ffccc6c2f
if save_data
	freq_save = filtered_omegas
	amp_save = filtered_amps
	phase_save = phases
	time_save = shifted_time
	if bath == "with bathymetry"
		jldsave("irregWaveData_withBathy.jld2";amp_save,freq_save,phase_save,time_save)
	else
		jldsave("irregWaveData_noBathy.jld2";amp_save,freq_save,phase_save,time_save)
	end
end

# ╔═╡ b0da79bb-7572-4ca1-a84e-d26a1c056b94
begin
	freqs = 0:0.1:10
	h = 0.3
	ks = computeWavenumber.(freqs,(h,))
	lambdas = 2*pi ./ks
	f3 = Figure()
	ax = Axis(f3[1,1], yticks=WilkinsonTicks(6,k_min=5))
	lines(freqs,lambdas)
	lines!(filtered_omegas,filtered_amps*30/0.002)
	current_figure()
end

# ╔═╡ dcdb49a6-707b-4d71-9cac-a8c0269eaff4
@bind phase_reg PlutoUI.Slider(0:0.01:2*pi,show_value=true)

# ╔═╡ 165e7b0b-0131-4446-8835-c5babf5c6c18
begin
	reg_wave = SimpleWave(0.005,2*pi*0.35,phase_reg,hasFadeIn=false)
	reg_eta = -1/9.81*[analyticPotential_dt(0,0,t,0.3,reg_wave) for t in shifted_time]
end

# ╔═╡ aa17f115-1e0f-4574-b819-ddd29cf3575d
begin
	savefigure=false
	fig5 = Figure(size=(500,300),font="CMU")
	ax5 = Axis(fig5[1,1], xlabel= L"$t$\,/s",ylabel=L"$\eta$\,/m")
	lines!(ax5,shifted_time,selected_eta[1:end-1],label="measurement",color=1, colormap=:viridis, colorrange=(1,10))
	lines!(ax5,shifted_time,reg_eta,label="regular wave",color=5, colormap=:viridis, colorrange=(1,10))
	lines!(ax5,shifted_time,irreg_eta,label="irregular wave",color=9, colormap=:viridis, colorrange=(1,10),linestyle=:dash)
	axislegend(position = :lt)
	if savefigure
		save("C:/Users/chris/Desktop/Masterarbeit/Text/masterthesis/clean-math-thesis/images/incomingWaves.svg",fig5)
	end
	fig5
end

# ╔═╡ Cell order:
# ╟─ac91d146-a4d4-44a7-ace7-fdb770f0b367
# ╠═7e8dab8d-fc6e-4a67-8ba1-87e48955b627
# ╟─2deecc02-48f9-46d3-8ac2-fba747562cbb
# ╟─1834c732-00e8-4cbf-9869-7853f9174af3
# ╟─58d9e3eb-2187-4eee-8ea1-457429b38fda
# ╟─cad908ae-1fa3-4a04-ae91-282d34ef4562
# ╟─902d2a0a-cd9c-4f56-a04e-3c155fcea3e4
# ╟─26daba59-0f1d-4dbe-8e70-1cc78ddb2f5e
# ╟─e4cc7114-3b89-48e7-847a-06bd961e10d8
# ╟─66c4193d-fbc2-4a67-84dd-ef18652171b6
# ╟─f8c7bc9c-76ab-4c4c-aff6-238e622a21bf
# ╟─b5634eba-2bc9-4de5-a13e-afb4a7cb38e7
# ╟─8ae9e601-bbdd-4e72-b8b0-13218d77e91c
# ╟─a0af295d-2850-41f6-a6bc-e5576b957bd7
# ╟─1ecbc240-6615-4249-8032-f569f1fb8524
# ╟─5f652052-e23b-409d-8776-074e13a49699
# ╟─0860d49f-33b1-4385-8aad-5d06c44715d4
# ╟─3c85ca95-00a7-4307-ae5f-9c6f9833d2f9
# ╠═d4cb9ea5-7b1e-47d6-ba81-d434c4624727
# ╟─3edc7820-70e5-4128-a472-a7ab7a1c095f
# ╠═4f9c2646-6552-4cfa-9ff8-b8eefba7de86
# ╠═3ac0bd7a-0248-4e18-b64e-32ec0d9c7404
# ╠═7529d603-1211-464c-903c-e45ffccc6c2f
# ╟─b0da79bb-7572-4ca1-a84e-d26a1c056b94
# ╠═4bc4bd94-17c3-4181-a70e-f02f9cad4f02
# ╠═dcdb49a6-707b-4d71-9cac-a8c0269eaff4
# ╠═165e7b0b-0131-4446-8835-c5babf5c6c18
# ╠═aa17f115-1e0f-4574-b819-ddd29cf3575d
