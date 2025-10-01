#This module has all the remeshing and adaptive derivatives
#uses proper fornberg algorithm to ge the derivative weights
#This writes the derivatives explicitly as loops - This is much faster
module fornberg
using  Dierckx #needed for splining

#a struct that stores the current state
mutable struct state
    h::Vector{Float64} #thickness
    c::Vector{Float64} #centerline
    coords::Vector{Float64} #coordinates
    Δx::Vector{Float64} #Vector of spacings
    t ::Float64 #times
end

#takes in coordinates and a list of data points, computes the halfwidth
function halfwidth(dat::Vector,coords::Vector,dxmax)
    min_index = argmin(dat) #min index
    min_dat = dat[min_index] #value of the minimum
    min_coord = coords[min_index] #coordinate of minimum
    width_index = findfirst(x->x<2*min_dat,dat) #the half width index may be found as the first point less than twice the min value
    
    #if we it hasn't decreased too much then return dxmax
    if isnothing(width_index)==true
        width=20*dxmax
        dxmin=dxmax
    else
        width = abs(coords[width_index]-min_coord) #half width
        dxmin = min(width/20,dxmax)
    end
    #if the width is too big
    if width> 20*dxmax 
        width=20*dxmax
        dxmin=dxmax
    end
    #return the minimum position, width, and new minimum spacing
    return min_coord, width, dxmin 
end



function make_grid(min_coord,dxmin,dxmax,L,ratio,width)
    #uniformly mesh central region around the singularity
    #we do it in these two steps to ensure we include a point at min_coord
    coords = collect(min_coord:dxmin:min_coord+4*width)
    prepend!(coords,collect(min_coord-4*width:dxmin:min_coord-dxmin))
    dx0=copy(dxmin)
    #slowly change the grid outside the central region and mesh the part around
    #mesh the region to the right of the pinch
    while coords[end]<L-2*dxmax
        #calculate the new spacing
        dx0 = min((coords[end]-coords[end-1])*ratio,dxmax)
        #add the new spacing
        push!(coords,coords[end]+dx0)
    end
    #mesh the left region
    while coords[1]>2*dxmax
        #calculate the new spacing
        dx0 = min((coords[2]-coords[1])*ratio,dxmax)
        #add the new spacing
        pushfirst!(coords,coords[1]-dx0)
    end
    #push!(coords, L-coords[1])
    prepend!(coords,[0.0,coords[1]/2])
    Δx = push!(diff(coords),mod(coords[1]-coords[end],L))
    return coords,Δx
end



#remeshes the old data onto a new grid using Dierckx splining
function spline(data::Vector,old_coords::Vector,coords::Vector,L::Float64)
    #Need to use a much higher order spline otherwise we accumulate too much error when taking high order
    #derivatives
    #spline = interpolate(old_coords,data,BSplineOrder(7);Periodic)
    #return spline.(coords)
    
    #In the spline formula we push points to the end
    #these are needed for periodicity as Dierckx assumed the end and start are the same
    c_data=copy(data)
    c_old = copy(old_coords)
    push!(c_data,data[1])
    push!(c_old,L)
    spline = Dierckx.Spline1D(c_old,c_data ;k=5, periodic=true) #ensure a periodic domain 
    return Dierckx.evaluate(spline,coords) #evaluate on the new coordinates
end

#takes in the current state and remeshes it
function remesh(state,dxmax,L,ratio)
    min_coord, width, dxmin  = halfwidth(state.h,state.coords,dxmax)
    coords, Δx = make_grid(min_coord,dxmin,dxmax,L,ratio,width)
    new_h = spline(state.h,state.coords,coords,L)
    new_c = spline(state.c,state.coords,coords,L)
    return coords,Δx,new_h,new_c
end

#function that computes the first derivative evaluated at x_i
function ∂1(f,Δx)
    n=length(f)
    [ (f[mod1(i+1,n)]-f[mod1(i-1,n)])/(Δx[i]+Δx[mod1(i-1,n)]) for i in 1:n]
end

#function that computes the first derivative evaluated at x_{i+1/2}
function ∂h1(f,Δx)
    n=length(f)
    [ (f[mod1(i+1,n)]-f[i])/(Δx[i]) for i in 1:n]
end

#function that computes the second derivative at x_i
function ∂2(f,Δx)
    n=length(f)
    [(2/(Δx[i]+Δx[mod1(i-1,n)]))*((f[mod1(i+1,n)]-f[i])/(Δx[i])-(f[i]-f[mod1(i-1,n)])/(Δx[mod1(i-1,n)])) for i in 1:n]
end



end



#=
using .fornberg,Plots, BenchmarkTools
L=2*π;
dxmax=5e-3;
uni_coords=collect(0:dxmax:L-dxmax);


thick(x)=2 .+ 1.9995*cos.(x);
thick3(x)=- 1.9995*sin.(x);
min_coord, width, dxmin = fornberg.halfwidth(thick.(uni_coords),uni_coords,dxmax)
coords,dxs=fornberg.make_grid(min_coord,dxmin,dxmax,L,1.02,width);
bla=fornberg.spline(thick.(uni_coords),uni_coords,coords,L);

plot(coords,bla)
plot(coords,fornberg.∂1((bla.^2).*fornberg.∂3(bla,dxs),dxs))
plot(coords,fornberg.∂3(bla,dxs))
plot(coords,fornberg.∂4(bla,dxs))

fornberg.∂4(bla,dxs)

@btime fornberg.spline(thick.(uni_coords),uni_coords,coords,L);


=#
