# redo with amp = 0.005

using JLD2, CairoMakie
using InhomogeneousBathymetryEuler
virColors=[1,9]
MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Regular.otf")
))
function linToCart(i_lin::Integer,N::Integer,M::Integer)
    A = zeros(N,M)
    return CartesianIndices(A)[i_lin].I[end:-1:1]
end
saveFigs = false

#data = load("examples/results/plottingScripts/relaxationVsDirectCompDataMoreShifts.jld2")
data = load("examples/results/plottingScripts/relaxationVsDirectCompData2.jld2")
etas = data["etas"]
energies = data["energies"]
wave = data["wave"]
lam = 2*pi/computeWavenumber(wave,0.3)
transs = data["transs"]
time_vec = data["time_vec"]
χss = data["χss"]
xss = [transs[1,1].x(χss[idx]) for idx = CartesianIndices(χss)]
gauss_shifts = data["gauss_shifts"]
Dχ = χss[1,1][2] - χss[1,1][1]
nShifts = size(etas,1)
nRx = size(etas,2)
gauss_shifts_str = [L"$ %$(round(shift,digits=1))$\,m" for shift in gauss_shifts]
gen_str = ["Direct Generation","L_RX = 1λ","L_RX = 2λ","L_RX = 4λ"]
peak_inds = [argmin(i -> abs(gauss_shifts[idx_s]-transs[idx_s,1].x(χss[idx_s,idx_w][i])), eachindex(χss[idx_s,idx_w])) for idx_s = 1:nShifts, idx_w = 1:nRx]
five_m_shift_ind = -241
#five_m_shift_ind = +100
behind_peak_inds = peak_inds .+ five_m_shift_ind
χs_phys = χss[1,1][230*2+1:end]
equi_inds_t = 1801:2001

eta_t_5behind_gauss = Matrix{Vector{Float64}}(undef,nShifts,nRx)
for idx_s = 1:nShifts
    for idx_w = 1:nRx
        eta_t_5behind_gauss[idx_s,idx_w] = [etas[idx_s,idx_w][t][behind_peak_inds[idx_s,idx_w]] for t=eachindex(time_vec)]
    end
end

figs_ = Figure[]
for (idx_fig,idx_s) = enumerate(1:nShifts)
    push!(figs_,plotSurfaceOverTime(etas[idx_s,1],time_vec,χss[idx_s,1],"direct"))
    #plotSurfaceOverTime!(figs_[idx_fig],etas[idx_s,2],time_vec,χss[idx_s,2],"rx 1")
    plotSurfaceOverTime!(figs_[idx_fig],etas[idx_s,3],time_vec,χss[idx_s,3],"rx 2")
    #plotSurfaceOverTime!(figs_[idx_fig],etas[idx_s,4],time_vec,χss[idx_s,4],"rx 4")
    ax = figs_[idx_fig].content[1]
    xlims!(ax,(-63,0))
    Label(figs_[idx_fig][begin-1,1:2],"shift = $(transs[idx_s,1].χ(gauss_shifts[idx_s]))")
end
#=
fig2 = Figure(size=(700,400))
ax21 = Axis(fig2[1,1],xlabel="t (s)", ylabel="η",title="Direct Wave Generation")
ax22 = Axis(fig2[1,2],xlabel="t (s)", ylabel="η",title="Relaxed Wave Generation, 1")
ax23 = Axis(fig2[2,1],xlabel="t (s)", ylabel="η",title="Relaxed Wave Generation, 2")
ax24 = Axis(fig2[2,2],xlabel="t (s)", ylabel="η",title="Relaxed Wave Generation, 4")
for idx_s = 2:nShifts
    lines!(ax21,time_vec[equi_inds_t],eta_t_5behind_gauss[idx_s,1][equi_inds_t],label="shift="*gauss_shifts_str[idx_s])
    lines!(ax22,time_vec[equi_inds_t],eta_t_5behind_gauss[idx_s,2][equi_inds_t],label="shift="*gauss_shifts_str[idx_s])
    lines!(ax23,time_vec[equi_inds_t],eta_t_5behind_gauss[idx_s,3][equi_inds_t],label="shift="*gauss_shifts_str[idx_s])
    lines!(ax24,time_vec[equi_inds_t],eta_t_5behind_gauss[idx_s,4][equi_inds_t],label="shift="*gauss_shifts_str[idx_s])
end
axislegend(ax21,position=:lt)
axislegend(ax22,position=:lt)
axislegend(ax23,position=:lt)
axislegend(ax24,position=:lt)
ylims!(ax21,-0.0075,0.0075)
ylims!(ax22,-0.0075,0.0075)
ylims!(ax23,-0.0075,0.0075)
ylims!(ax24,-0.0075,0.0075)
fig2 =#

max_eta_t_5behind_gauss = [maximum(abs.(eta_t_5behind_gauss[idx][equi_inds_t])) for idx = CartesianIndices(eta_t_5behind_gauss)]
#= curr_max_eta_t_5behind_gauss = Matrix{Vector{Float64}}(undef,nShifts,nRx)
for idx_s = 1:nShifts
    for idx_w = 1:nRx
        M=0
        curr_max_eta_t_5behind_gauss[idx_s,idx_w] = zeros(length(time_vec))
        for idx_t = 1:length(time_vec)
            M = abs(eta_t_5behind_gauss[idx_s,idx_w][idx_t]) > M ? abs(eta_t_5behind_gauss[idx_s,idx_w][idx_t]) : M
            curr_max_eta_t_5behind_gauss[idx_s,idx_w][idx_t] = M
        end
    end
end
figs2_ = Figure[]
for (idx_fig,idx_s) = enumerate(1:nShifts)
    push!(figs2_,Figure())
    ax = Axis(figs2_[idx_fig][1,1],xlabel="t (s)")
    for idx_w = [1,3]
        lines!(ax,time_vec,curr_max_eta_t_5behind_gauss[idx_s,idx_w],label=gen_str[idx_w])
    end
    axislegend(ax,position=:rb)
    ylims!(ax,0,0.01)
end =#


fig3 = Figure(size=(600,250))
ax3 = Axis(fig3[1,1:2],xlabel=L"$x_G$\,/m",ylabel=L"Amplitude\,/m$$",xticks=WilkinsonTicks(6,k_min=4))
label_ad = ["Direct Generation",nothing,"Relaxed Generation"]
for (i,idx_w) = enumerate([1,3])
    lines!(ax3,gauss_shifts,max_eta_t_5behind_gauss[:,idx_w],label=label_ad[idx_w],color=virColors[i], colormap=:viridis, colorrange=(1,10))
end
Legend(fig3[1,3],ax3)
#axislegend(ax3,position=:rb)
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\dirRXamps.svg",fig3)
end
fig3

v_p = wave.freq/wave.waveNumberAtInflow
t_inds_gauss_hit_by_wave = [argmin(t_ind->abs((5+gauss_shifts[idx_s])/v_p - time_vec[t_ind]),eachindex(time_vec)) for idx_s=1:nShifts]
fig4 = Figure(size=(600,450))
how_much_time = 800
axs4 = Axis[]
virColors=[1,5,9]
label_ad = ["Direct Wave Generation",L"Relaxed Wave Generation, $L_{RX}=1\lambda$",L"Relaxed Wave Generation, $L_{RX}=2\lambda$",L"Relaxed Wave Generation, $L_{RX}=4\lambda$"]
for idx_w = 1:4
    push!(axs4,Axis(fig4[linToCart(idx_w,2,2)...],xlabel=L"$t$\,/s",ylabel=L"$\hat{\eta}(x_G,t)$\,/m",yticks=[-0.005,0,0.005],title=label_ad[idx_w]))
    for (i,idx_s)=enumerate([1,4,8])
        t_inds_s = t_inds_gauss_hit_by_wave[idx_s]:t_inds_gauss_hit_by_wave[idx_s]+how_much_time
        lines!(axs4[idx_w],time_vec[t_inds_s].-time_vec[t_inds_gauss_hit_by_wave[idx_s]],eta_t_5behind_gauss[idx_s,idx_w][t_inds_s],label=gauss_shifts_str[idx_s],color=virColors[i], colormap=:viridis, colorrange=(1,10))
    end
    ylims!(axs4[idx_w],-0.0075,0.0075)
end
Legend(fig4[3,1:2],axs4[1],"Bathymetry Position",orientation=:horizontal)
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\dirRXtimeLocal.svg",fig4)
end
fig4

virColors=[1,9]
fig5 = Figure(size=(700,300))
ax5 = Axis(fig5[1,1],xlabel=L"$t$\,/s",ylabel=L"$\eta$\,/m",xticks=WilkinsonTicks(6,k_min=5))
lines!(ax5,time_vec[1:2001],eta_t_5behind_gauss[end-2,1][1:2001],label="Direct Generation",color=virColors[1], colormap=:viridis, colorrange=(1,10))
lines!(ax5,time_vec[1:2001],eta_t_5behind_gauss[end-1,3][1:2001],label="Relaxed Generation",color=virColors[2], colormap=:viridis, colorrange=(1,10))
axislegend(ax5,position=:rb)
if saveFigs
    save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\dirRXoverTime.svg",fig5)
end
fig5