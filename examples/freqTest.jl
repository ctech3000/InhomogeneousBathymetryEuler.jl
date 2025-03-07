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

# ╔═╡ 26daba59-0f1d-4dbe-8e70-1cc78ddb2f5e
begin
	eta = zeros(Float64,10000)
	time_sensor = zeros(Float64,10000)
	for heat_idx = 1:20
		filename = "C:\\Users\\chris\\.julia\\dev\\InhomogeneousBathymetryEuler.jl\\examples\\Data_Sensors\\Without_Bathymetry\\Heat$(heat_idx).txt"
		sensor = 1
		sensor_pos = [0.0,2.0,4.0,6.0]
		x = sensor_pos[sensor]
		z = 0.0
		global eta += convert(Vector{Float64},readdlm(filename)[2:end,sensor])/100
		global time_sensor = convert(Vector{Float64},readdlm(filename)[2:end,5])
	end
	eta /= 20
end

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

# ╔═╡ 1ecbc240-6615-4249-8032-f569f1fb8524
begin 
	lines(selected_time,selected_eta)
end

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
Select maximum k (number of components).
"""

# ╔═╡ 3c85ca95-00a7-4307-ae5f-9c6f9833d2f9
@bind max_k PlutoUI.Slider(eachindex(coeff), default = 60, show_value=true)

# ╔═╡ d4cb9ea5-7b1e-47d6-ba81-d434c4624727
begin
	filtered_amps = amps[1:max_k]
	filtered_omegas = omegas[1:max_k]
	filtered_cos_reco = [sum(filtered_amps.*cos.(filtered_omegas*t)) for t in shifted_time]
end;

# ╔═╡ 3edc7820-70e5-4128-a472-a7ab7a1c095f
begin
	fig3 = Figure(size = (600,900))
	ax3 = Axis(fig3[1,1],xlabel="t",title="reconstruction")
	lines!(ax3,shifted_time,selected_eta[1:end-1],label="orig.")
	lines!(ax3,shifted_time,filtered_cos_reco, label="filtered cos_reco")
	axislegend(position = :lb)
	ax3_5 = Axis(fig3[2,1],title="coefficients")
	lines!(ax3_5,filtered_omegas,filtered_amps)
	fig3
end

# ╔═╡ 4f9c2646-6552-4cfa-9ff8-b8eefba7de86
begin
	phases = zeros(size(filtered_amps))
	wave = IrregWave(filtered_amps[2:end],filtered_omegas[2:end],phases[2:end],hasFadeIn=false)
	comp_eta = -1/9.81*[analyticPotential_dt(0,0,t,0.3,wave) for t in shifted_time]
	fig4 = Figure()
	ax4 = Axis(fig4[1,1])
	lines!(ax4,shifted_time,selected_eta[1:end-1],label="orig.")
	lines!(ax4,shifted_time,filtered_cos_reco, label="filtered cos_reco")
	lines!(ax4,shifted_time,comp_eta.+amps[1],label="comp_eta")
	axislegend(position = :lt)
	fig4
end

# ╔═╡ 3ac0bd7a-0248-4e18-b64e-32ec0d9c7404
save_data = true

# ╔═╡ 7529d603-1211-464c-903c-e45ffccc6c2f
if save_data
	freq_save = filtered_omegas
	amp_save = filtered_amps
	phase_save = phases
	time_save = shifted_time
	jldsave("sensorFreqsWith0.jld2";amp_save,freq_save,phase_save,time_save)
end

# ╔═╡ Cell order:
# ╠═ac91d146-a4d4-44a7-ace7-fdb770f0b367
# ╠═7e8dab8d-fc6e-4a67-8ba1-87e48955b627
# ╠═2deecc02-48f9-46d3-8ac2-fba747562cbb
# ╠═26daba59-0f1d-4dbe-8e70-1cc78ddb2f5e
# ╟─e4cc7114-3b89-48e7-847a-06bd961e10d8
# ╠═66c4193d-fbc2-4a67-84dd-ef18652171b6
# ╟─f8c7bc9c-76ab-4c4c-aff6-238e622a21bf
# ╟─b5634eba-2bc9-4de5-a13e-afb4a7cb38e7
# ╠═8ae9e601-bbdd-4e72-b8b0-13218d77e91c
# ╟─1ecbc240-6615-4249-8032-f569f1fb8524
# ╟─5f652052-e23b-409d-8776-074e13a49699
# ╟─0860d49f-33b1-4385-8aad-5d06c44715d4
# ╟─3c85ca95-00a7-4307-ae5f-9c6f9833d2f9
# ╠═d4cb9ea5-7b1e-47d6-ba81-d434c4624727
# ╟─3edc7820-70e5-4128-a472-a7ab7a1c095f
# ╠═4f9c2646-6552-4cfa-9ff8-b8eefba7de86
# ╠═3ac0bd7a-0248-4e18-b64e-32ec0d9c7404
# ╠═7529d603-1211-464c-903c-e45ffccc6c2f
