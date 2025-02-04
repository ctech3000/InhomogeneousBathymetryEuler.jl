using GLMakie

f = Figure()
ax = Axis3(f[1, 1], xlabel = "χ", ylabel = "σ", zlabel="ϕ")
surface!(ax,χs,σs,all_phis[200])
f

f2 = Figure()
ax2 = Axis(f2[1, 1], xlabel="χ", ylabel="η")
lines!(ax2,χs,all_etas[80])
f2

GLMakie.activate!()
amp = wave.amp
f3 = Figure()
#ax3 = Axis(f3[1,1],xlabel="χ",ylabel="η")
ax3_x_D = Axis(f3[1,1], yticks=[-1.5*amp,-amp,amp,1.5*amp], xticks=([trans.χ(domain.x_D)],["χ_D"]),xlabel="χ", ylabel="η")
time_ind = 1:length(time_vec)
t_sl = Slider(f3[2, 1], range = time_ind, startvalue = 1)
eta = lift(t_sl.value) do t
    all_etas[t]
end
lines!(χs,eta)
#ylims!(ax3,-0.1,0.1)
ylims!(ax3_x_D,-1.5*amp,1.5*amp)
f3

f4 = Figure()
phi_amp = g*amp/wave.freq
#ax4 = Axis(f4[1,1],xlabel="χ",ylabel="η")
ax4_x_D = Axis(f4[1,1], yticks=[-1.5*phi_amp,-phi_amp,phi_amp,1.5*phi_amp], xticks=([trans.χ(domain.x_D)],["χ_D"]),xlabel="χ", ylabel="ϕ")
time_ind = 1:length(time_vec)
t_sl = Slider(f4[2, 1], range = time_ind, startvalue = 1)
eta = lift(t_sl.value) do t
    all_phis[t][:,1]
end
lines!(χs,eta)
#ylims!(ax4,-0.1,0.1)
ylims!(ax4_x_D,-1.5*phi_amp,1.5*phi_amp)
f4


#= t = Observable(1)
f3 = Figure()
data = @lift(all_etas[$t])
ax3 = Axis(f3[1,1],xlabel="χ",ylabel="η")
ls = labelslider!(f3, "t", time_ind)
f[2,1] = ls.layout
connect!(t, ls.slider.value)
f3 =#

