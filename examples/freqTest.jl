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
using InhomogeneousBathymetryEuler, PlutoUI, GLMakie, FFTW, DelimitedFiles

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
	filename = "C:\\Users\\chris\\.julia\\dev\\InhomogeneousBathymetryEuler.jl\\examples\\Data_Sensors\\Without_Bathymetry\\Heat1.txt"
	sensor = 1
	eta = convert(Vector{Float64},readdlm(filename)[2:end,sensor])/100
	time = convert(Vector{Float64},readdlm(filename)[2:end,5])
end

# ╔═╡ e4cc7114-3b89-48e7-847a-06bd961e10d8
md"""
Select start time.
"""

# ╔═╡ 66c4193d-fbc2-4a67-84dd-ef18652171b6
begin
	@bind t_start PlutoUI.Slider(time, default = time[1], show_value=true)
end

# ╔═╡ f8c7bc9c-76ab-4c4c-aff6-238e622a21bf
md"""
Select end time.
"""

# ╔═╡ b5634eba-2bc9-4de5-a13e-afb4a7cb38e7
@bind t_end PlutoUI.Slider(time, default = length(time), show_value=true)

# ╔═╡ 8ae9e601-bbdd-4e72-b8b0-13218d77e91c
begin
	t_start_ind = findall(x->x==t_start,time)[1]
	t_end_ind = findall(x->x==t_end,time)[1]
	selected_time_inds = t_start_ind:t_end_ind
	selected_time = time[selected_time_inds]
	selected_eta = eta[selected_time_inds]
end;

# ╔═╡ 1ecbc240-6615-4249-8032-f569f1fb8524
begin 
	lines(selected_time,selected_eta)
end

# ╔═╡ fa87a573-6ff9-4a90-9c96-f1fa93356308
begin 
	Fs = 100
	N = length(selected_time_inds)
	freqs = 0:Fs/N:(N-1)*Fs/N
	coeff = (fft(selected_eta))
end;

# ╔═╡ cd2f1498-7cd9-467e-a094-a319c8ac3d49
md"""
Select start frequency.
"""

# ╔═╡ f845464e-5c48-479b-86c6-f90115692986
@bind start_frequency PlutoUI.Slider(freqs,default = freqs[1], show_value = true)

# ╔═╡ f8b9004d-22ae-49b5-a77a-3893577448d3
md"""
Select end frequency.
"""

# ╔═╡ cd663bf7-07ca-4bb5-aff9-8f832b8c7215
@bind end_frequency PlutoUI.Slider(freqs,default = freqs[end], show_value = true)

# ╔═╡ 672f4203-9540-405b-8834-e381fdf207e4
begin
	freq_start_ind = findall(x->x==start_frequency,freqs)[1]
	freq_end_ind = findall(x->x==end_frequency,freqs)[1]
	selected_freq_inds = freq_start_ind:freq_end_ind
	selected_freqs = freqs[selected_freq_inds]
	selected_coeffs = coeff[selected_freq_inds]
end;

# ╔═╡ aa7d7d32-8f5c-4978-b070-ee1ac080f564
begin
	f1 = Figure()
	ax1 = Axis(f1[1,1],title="selected amps")
	lines!(ax1,2*pi*selected_freqs,real.(selected_coeffs))
	f1
end

# ╔═╡ c3024e1c-f7cf-47ea-844a-2f6936b3ca80
begin
	filtered_coeff = zeros(eltype(coeff),size(coeff))
	filtered_coeff[selected_freq_inds] = coeff[selected_freq_inds]
	filtered_coeff[end.-(selected_freq_inds.-1)] = coeff[end.-(selected_freq_inds.-1)]
	reconstruction = real.(ifft(filtered_coeff))
	f = Figure()
	ax = Axis(f[1,1],title="Fourier Representation")
	lines!(ax,selected_time,selected_eta,label="original")
	lines!(ax,selected_time,reconstruction,label="filtered")
	axislegend(position = :lb)
	f
end

# ╔═╡ 0211f9c6-8d8d-4b90-afc7-ed123b32cc56
begin
	wave_amps = real.(vcat(coeff[selected_freq_inds], coeff[end.-(selected_freq_inds.-1)]))
	wave_freqs = 2*pi*vcat(freqs[selected_freq_inds],freqs[end.-(selected_freq_inds.-1)])
	lines(wave_freqs,wave_amps)
end

# ╔═╡ Cell order:
# ╠═ac91d146-a4d4-44a7-ace7-fdb770f0b367
# ╠═7e8dab8d-fc6e-4a67-8ba1-87e48955b627
# ╠═2deecc02-48f9-46d3-8ac2-fba747562cbb
# ╟─26daba59-0f1d-4dbe-8e70-1cc78ddb2f5e
# ╟─e4cc7114-3b89-48e7-847a-06bd961e10d8
# ╟─66c4193d-fbc2-4a67-84dd-ef18652171b6
# ╟─f8c7bc9c-76ab-4c4c-aff6-238e622a21bf
# ╟─b5634eba-2bc9-4de5-a13e-afb4a7cb38e7
# ╟─8ae9e601-bbdd-4e72-b8b0-13218d77e91c
# ╟─1ecbc240-6615-4249-8032-f569f1fb8524
# ╠═fa87a573-6ff9-4a90-9c96-f1fa93356308
# ╟─cd2f1498-7cd9-467e-a094-a319c8ac3d49
# ╟─f845464e-5c48-479b-86c6-f90115692986
# ╟─f8b9004d-22ae-49b5-a77a-3893577448d3
# ╟─cd663bf7-07ca-4bb5-aff9-8f832b8c7215
# ╟─672f4203-9540-405b-8834-e381fdf207e4
# ╠═aa7d7d32-8f5c-4978-b070-ee1ac080f564
# ╠═c3024e1c-f7cf-47ea-844a-2f6936b3ca80
# ╠═0211f9c6-8d8d-4b90-afc7-ed123b32cc56
