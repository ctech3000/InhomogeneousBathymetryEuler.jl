using CairoMakie
using InhomogeneousBathymetryEuler
CairoMakie.activate!()
vifColors=[1,5,9]

bathPoints = collect(1.8:0.001:3.2)
bath = Bathymetry(bathPoints,"Gauss",shift=2.5)

fig = Figure(size=(400,200))
ax = Axis(fig[1,1],xlabel=L"$x$\,/m",ylabel=L"$b$\,/m")
lines!(ax,bathPoints,eval_bath(bath,bathPoints),color=1, colormap=:viridis, colorrange=(1,10))
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\bathymetry.svg",fig)
fig

bathPoints2 = collect(0:0.001:1.5)
baths = Bathymetry[]
push!(baths,Bathymetry(bathPoints2,"Gauss",shift=0.75,bHeight=0.375,depth=-1.5))
push!(baths,Bathymetry(bathPoints2,"Gauss",shift=0.75,bHeight=0.75,depth=-1.5))
push!(baths,Bathymetry(bathPoints2,"Gauss",shift=0.75,bHeight=1.125,depth=-1.5))
bathNames = [L"\hat{b}_1",L"\hat{b}_2",L"\hat{b}_3"]

fig2 = Figure(size=(400,250))
ax2 = Axis(fig2[1,1],xlabel=L"$x$\,/m",ylabel=L"$b$\,/m")
for idx_b = 1:3
    lines!(ax2,bathPoints2,eval_bath(baths[idx_b],bathPoints2),color=virColors[idx_b], colormap=:viridis, colorrange=(1,10),label=bathNames[idx_b])
end
axislegend(ax2,position=:lt)
save("C:\\Users\\chris\\Desktop\\Masterarbeit\\Text\\masterthesis\\clean-math-thesis\\images\\bathymetryManSol.svg",fig2)
fig2
