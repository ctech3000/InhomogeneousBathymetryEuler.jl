### A Pluto.jl notebook ###
# v0.20.0

using Markdown
using InteractiveUtils

# ╔═╡ ccdcedb0-1e16-11f0-27d4-fb65f20d759e
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate(Base.current_project())
    using Revise # Revise must be in your LOAD_PATH
end;

# ╔═╡ 287e4ed7-c53c-4d31-a6e8-af9f085e96e8
using InhomogeneousBathymetryEuler, PlutoUI, GLMakie, JLD2;

# ╔═╡ fdc470cb-5dc4-46ae-8d24-861f9d094139
html"""
<style>
input[type*="range"] {
	width: 100%;
}
</style>
"""

# ╔═╡ 8910b1db-d393-4d16-bfde-3cbc57f5c9d3
begin
	d_sine_flat = load("presentationData/data_sine_flat.jld2")
	eta_sine_flat = d_sine_flat["all_etas"]
	bath_flat = d_sine_flat["bath_val"]
	time_vec_sine = d_sine_flat["time_vec"]
	sensors_sine_flat = d_sine_flat["sensors"]
	xs_sine = d_sine_flat["xs"]
	d_sine_gauss = load("presentationData/data_sine_gauss.jld2")
	eta_sine_gauss = d_sine_gauss["all_etas"]
	bath_gauss = d_sine_gauss["bath_val"]
	bathPoints = d_sine_gauss["bathPoints"]
	sensors_sine_gauss = d_sine_gauss["sensors"]
	phys_inds_sine = 801:2401
	sampled_phys_inds = 801:8:2401
	trans_bath_gauss = (bath_gauss[1:7:2668].+0.3)*0.005/0.2.-0.015
	bathPoints_sampled = bathPoints[1:7:2668]
end;

# ╔═╡ 90c4a3d1-6e4b-42e0-84eb-03928b10d964
md"""
## Sine wave over flat bottom
"""

# ╔═╡ 275296cd-22bb-4895-a2e9-5f86b7e00f1a
begin
	time1 = Observable(1)

	ys1 = @lift(eta_sine_flat[$time1][sampled_phys_inds])
	
	fig1 = Figure()
	axis1 = Axis(fig1[1,1],xlabel="x/m",ylabel="η/m",title = @lift("t = $(round(time_vec_sine[$time1], digits = 1))s"))
	lines!(axis1,xs_sine[sampled_phys_inds], ys1, color = :blue, linewidth = 2.5)
	ylims!(axis1,(-0.006,0.006))
	#lines!(axis1,xs_sine[sampled_phys_inds],-0.015*ones(size(sampled_phys_inds)),color=:brown, linewidth = 1.5)
	
	    
	timestamps = range(1, 1001, step=1)
	
	GLMakie.Makie.Record(fig1, timestamps;
	        framerate = 80) do t
	        time1[] = t
	end
end

# ╔═╡ ab549aab-1935-484d-9f06-09f00edcd9a6
md"""
# Sine wave over gaussian bathymetry
"""

# ╔═╡ d1b0d20a-d473-4aae-bb8d-c152e42d6da4
begin
	time2 = Observable(1)

	ys2 = @lift(eta_sine_gauss[$time2][sampled_phys_inds])
	
	fig2 = Figure()
	axis2 = Axis(fig2[1,1],xlabel="x/m",ylabel="η/m",title = @lift("t = $(round(time_vec_sine[$time2], digits = 1))s"),yticks=[-0.01,-0.005,0,0.005,0.01])
	lines!(axis2,xs_sine[sampled_phys_inds], ys2, color = :blue, linewidth = 2.5)
	ylims!(axis2,(-0.016,0.01))
	lines!(axis2,bathPoints_sampled,trans_bath_gauss,color=:brown, linewidth = 1.5)
	timestamps2 = range(1, 2001, step=1)
	GLMakie.Makie.Record(fig2, timestamps2;
	        framerate = 40) do t
	        time2[] = t
	end
end

# ╔═╡ 30db3ddc-c4a2-4026-a1a5-f955c063f481
md"""
## Difference in surface displacement
"""

# ╔═╡ 81788ac6-0d0c-4f0b-aa30-98999a52f0b5
begin
	time4 = Observable(1)
	sampled_phys_inds_rx = 801:8:3201
	ys4 = @lift(eta_sine_gauss[$time4][sampled_phys_inds_rx].-eta_sine_flat[$time4][sampled_phys_inds_rx])
	
	fig4 = Figure()
	axis41 = Axis(fig4[1,1],xlabel="x/m",ylabel="η_d/m",title = @lift("t = $(round(time_vec_sine[$time4], digits = 1))s"))
	lines!(axis41,xs_sine[sampled_phys_inds_rx], ys4, color = :blue, linewidth = 2.5)
	ylims!(axis41,(-0.01,0.01))
	#lines!(axis41,bathPoints_sampled,trans_bath_gauss,color=:brown, linewidth = 1.5)
	
	timestamps4 = range(281, 2001, step=1)
	GLMakie.Makie.Record(fig4, timestamps4;
	        framerate = 40) do t
	        time4[] = t
	end
end

# ╔═╡ ca056942-167c-494f-847a-9af0e3c5020c
md"""
## Surface displacement over time with and without bathymetry
"""

# ╔═╡ daf7448a-986c-422e-80d4-bb5dcf1cdee5
begin
	fig3 = Figure(size=(800,800))
	ax31 = Axis(fig3[1,1],xlabel = "t/s",ylabel = "η/m",title="η at x = 8m")
	lines!(ax31,time_vec_sine,sensors_sine_flat.data[2],label="no bath.",linewidth = 2.5)
	lines!(ax31,time_vec_sine,sensors_sine_gauss.data[2],label="with bath.",linewidth = 2.5, color = GLMakie.Makie.wong_colors()[6])
	xlims!(7,42)

	ax32 = Axis(fig3[2,1],xlabel = "t/s",ylabel = "η/m",title="η at x = 12m")
	lines!(ax32,time_vec_sine,sensors_sine_flat.data[3],label="no bath.",linewidth = 2.5)
	lines!(ax32,time_vec_sine,sensors_sine_gauss.data[3],label="with bath.",linewidth = 2.5,color = GLMakie.Makie.wong_colors()[6])
	Legend(fig3[3,1],fig3.content[1],orientation=:horizontal)
	xlims!(7,42)
	fig3
end

# ╔═╡ 0a5f8177-5731-440b-8eba-dfc72af10cf3
begin
	d_irreg = load("presentationData/data_irreg.jld2")
	eta_irreg = d_irreg["all_etas"]
	xs_irreg = d_irreg["xs"]
	d_irreg_2 = load("presentationData/data_irreg_2.jld2")
	eta_irreg_2 = d_irreg_2["all_etas"]
	phys_inds_irreg = 801:1601
end;

# ╔═╡ 93daf9ba-b8b3-475c-8513-00485f7cc3a8
# ╠═╡ disabled = true
#=╠═╡
begin
	time5 = Observable(1)

	ys5 = @lift(eta_irreg[$time5][phys_inds_irreg])
	ys5_2 = @lift(eta_irreg_2[$time5][phys_inds_irreg])
	
	fig5 = Figure()
	axis5 = Axis(fig5[1,1],xlabel="x/m",ylabel="η/m",title = @lift("t = $(round(time_vec_sine[$time5], digits = 1))s"))
	lines!(axis5,xs_irreg[phys_inds_irreg], ys5, color = :blue, linewidth = 2.5)
	lines!(axis5,xs_irreg[phys_inds_irreg], ys5_2, color = :green, linewidth = 2.5)
	ylims!(axis5,(-0.01,0.01))	
	    
	timestamps5 = range(1, 2001, step=1)
	
	GLMakie.Makie.Record(fig5, timestamps5;
	        framerate = 80) do t
	        time5[] = t
	end
end
  ╠═╡ =#

# ╔═╡ 8baaa107-3502-46d3-94fa-6afbfb80565c
begin
	data = load("presentationData/dataComp_flat.jld2")
	time_vec = data["time_vec"]
	eta_flat = data["all_etas"]
	data = load("presentationData/dataComp_gauss.jld2")
	eta_gauss = data["all_etas"]
	xs = 25:-0.0125:-20
	time10 = Observable(1)
	ys10 = @lift(eta_flat[$time10].-eta_gauss[$time10])
	
	fig10 = Figure()
	axis101 = Axis(fig10[1,1],xlabel="x/m",ylabel="η_d/m",title = @lift("t = $(round(time_vec_sine[$time10], digits = 1))s"))
	lines!(axis101,xs, ys10, color = :blue, linewidth = 2.5)
	ylims!(axis101,(-0.01,0.01))
	#lines!(axis41,bathPoints_sampled,trans_bath_gauss,color=:brown, linewidth = 1.5)
	
	timestamps10 = range(1, length(time_vec), step=1)
	GLMakie.Makie.Record(fig10, timestamps10,framerate = 40) do t
	        time10[] = t
	end
end

# ╔═╡ Cell order:
# ╟─ccdcedb0-1e16-11f0-27d4-fb65f20d759e
# ╟─287e4ed7-c53c-4d31-a6e8-af9f085e96e8
# ╟─fdc470cb-5dc4-46ae-8d24-861f9d094139
# ╟─8910b1db-d393-4d16-bfde-3cbc57f5c9d3
# ╟─90c4a3d1-6e4b-42e0-84eb-03928b10d964
# ╟─275296cd-22bb-4895-a2e9-5f86b7e00f1a
# ╟─ab549aab-1935-484d-9f06-09f00edcd9a6
# ╟─d1b0d20a-d473-4aae-bb8d-c152e42d6da4
# ╟─30db3ddc-c4a2-4026-a1a5-f955c063f481
# ╠═81788ac6-0d0c-4f0b-aa30-98999a52f0b5
# ╟─ca056942-167c-494f-847a-9af0e3c5020c
# ╟─daf7448a-986c-422e-80d4-bb5dcf1cdee5
# ╟─0a5f8177-5731-440b-8eba-dfc72af10cf3
# ╠═93daf9ba-b8b3-475c-8513-00485f7cc3a8
# ╠═8baaa107-3502-46d3-94fa-6afbfb80565c
