cd("/Users/lukeneville/Documents/current/chiral_sims/fornberg_v5/");
include("fornberg.jl")
include("chiral_fornberg.jl")


using .fornberg, .chiral_forn, Plots, LaTeXStrings, JLD2,CurveFit,CairoMakie,LsqFit;
#Plot defaults
default(framestyle=:box,grid=false,
guidefont=font(12,"Computer Modern"),tickfont=font(12,"Computer Modern"),legendfont=font(12,"Computer Modern"),fg_legend = :transparent,dpi=300,
lw=2,legend=:false,linecolor=:black);
xtic = 10.0 .^ collect(range(-12,2,step=1));

dxmax = 1e-2; #maximum grid spacing
dtmax=1e-2;
L=Float64(π); #length of the domain
res=1e-6; #error requirement
ratio=1.05; #ratio between grid spacings
old_coords = collect(0:dxmax:L-dxmax); #uniform coords. Go one from dx at the end to avoid double counting in PBCs
Δxs = push!(diff(old_coords),L-old_coords[end]);
thickness = .5 .+  .1*cos.(old_coords*2*π/L);
centerline=.1*sin.(old_coords*2*π/L);



# initialise the state and data file
state0=fornberg.state(thickness,centerline,old_coords,Δxs,0.);
data=[[state0,dtmax]]; 

#@save "data_2.jld2" data



dtmin=1e-5;
chiral_forn.run(data,res,4,dxmax,L,ratio,dtmax,dtmin);



data[end][1].t,data[end][2]
minimum(data[end][1].h)
fornberg.halfwidth(data[end][1].h,data[end][1].coords,dxmax)


Plots.plot(data[1][1].coords,data[1][1].c.+data[1][1].h,linecolor=:red);
Plots.plot!(data[1][1].coords,data[1][1].c.-data[1][1].h,linecolor=:red);
Plots.plot!(data[end][1].coords,data[end][1].c.+data[end][1].h);
Plots.plot!(data[end][1].coords,data[end][1].c.-data[end][1].h)


Plots.plot(data[1][1].coords,data[1][1].c + data[1][1].h,ylims=(-5,5),xlabel=L"x",linecolor=:lightgray);
Plots.plot!(data[1][1].coords,data[1][1].c - data[1][1].h,linecolor=:lightgray)
Plots.plot!(data[end-2][1].coords,data[end-2][1].c - data[end-2][1].h, fillrange = data[end-2][1].c  .+ data[end-2][1].h, fillalpha = 0.7,label=:false,c=:darkorange2,linecolor=:darkorange3)
Plots.plot!(data[end-1][1].coords,data[end-1][1].c + data[end-1][1].h,linecolor=:darkorange3)
#savefig("pinch.pdf")




#extract important things
times=[data[i][1].t for i in 1:length(data)];
dc=[fornberg.∂1(data[i][1].c,data[i][1].Δx) for i in 1:length(data)];
thicks=[minimum(data[i][1].h) for i in 1:length(data)];
widths=[fornberg.halfwidth(data[i][1].h,data[i][1].coords,1e-2)[2] for i in 1:length(data)];
tsing=times[end]+7e-6;
tp=tsing.-times;
αt= fornberg.∂1(log.(thicks)[1:end-1],diff(log.(tp)));
Plots.plot(tp[1:end],thicks[1:end],yaxis=:log,xaxis=:log,xticks=xtic,yticks=xtic)


minimum(dc[end])


index=findfirst(x->x<2e-2,tp);
tfit = power_fit(tp[index:end-1],thicks[index:end-1])
Plots.plot(tp[1:end-1],thicks[1:end-1],yaxis=:log,xaxis=:log,xticks=xtic,ylabel="min(h)",xlabel=L"t'=t_0-t");
Plots.plot!(tp[index:end-1],tfit[1]*tp[index:end-1].^tfit[2],linecolor=:red,linestyle=:dash,yticks=xtic,minorticks=10)
#savefig("thickness_scaling.pdf")


wfit=power_fit(tp[index:end-1],widths[index:end-1])
Plots.plot(tp[1:end-1],widths[1:end-1],yaxis=:log,xaxis=:log,xticks=xtic,minorticks=10,ylabel="max(Δx)",xlabel=L"t'=t_0-t")
Plots.plot!(tp[index:end-1],1*wfit[1]*tp[index:end-1].^wfit[2],linecolor=:red,linestyle=:dash)


tfit[2]


(1+1.25)/4