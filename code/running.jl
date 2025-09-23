include("fornberg.jl")
include("chiral_fornberg.jl")

#need to have everything in the same folder

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

#initial thickness and centerline
thickness = .5 .+  .1*cos.(old_coords*2*π/L); 
centerline=.1*sin.(old_coords*2*π/L);



# initialise the state and data file
state0=fornberg.state(thickness,centerline,old_coords,Δxs,0.);
data=[[state0,dtmax]]; 




#below runs until the min time step is 1e-5 and the total time is 4.
#These can also be changed after run. I.e. if you run for a total time for 4 but want to keep going to 5, just change the number 
#and evaluate this line again
dtmin=1e-5; 
chiral_forn.run(data,res,4,dxmax,L,ratio,dtmax,dtmin);

#this saves the data to a JLD2 file. 
#@save "data.jld2" data




tfit[2]


(1+1.25)/4
