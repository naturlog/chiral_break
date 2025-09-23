module chiral_forn #simulates the chiral equations with the proper fornberg derivatives
include("fornberg.jl")
using Plots, .fornberg, NonlinearSolve,LinearSolve,SparseConnectivityTracer, ADTypes

#Computes the rhs of the chiral eqns
#combined into one vector, with c on top of h
function rhs(h,c,Δx)
    #mat is the vector of derivative matrices
    #compute first and second derivatives
    dc = fornberg.∂1(c,Δx)
    dh =fornberg.∂1(h,Δx)
    d2c = fornberg.∂2(c,Δx)
    d2h =  fornberg.∂2(h,Δx)
    #compute the rhs terms
    rhsc = (d2c./(h)) .- (4*dh./h)
    rhsh = fornberg.∂2(-4*h.*dc .+ (dc.^2 .+ dh.^2 .-2*h.*d2h)./2, Δx)
    #rhsh = fornberg.∂2(-2*h.*dc .+ dh.^2 .-2*h.*d2h, Δx) .+2*fornberg.∂1(dc.*d2c,Δx)
    #combine the terms into one big vector
    return vcat(rhsc,rhsh)
end

#compute an initial guess using a forward euler step
#return a vector that combines both
function guess(state,dt)
    vcat(state.c, state.h) .+ dt*rhs(state.h,state.c,state.Δx)
end

#this function (du) is zero to solve the problem
#according with the notation in the Nonlinearsolve package, u is what we want to find, and
#p are the parameters/ numbers that we use for it
function zero_func(du,u,p)
    n = div(length(u),2)
    #split the u vector into centerline and thickness
    c = u[1:n] 
    h = u[n+1:end]
    dt = p[2] #we always stick dt as the second parameter

    #p[1] is the old solution
    #p[3] is the vector of derivative matrices

    du .= u .- p[1] - dt*rhs(h,c,p[3])
end

#We use a fully implicit solving method
#solve using Newton's method and the initial guess computed earlier
#according with the notation in the Nonlinearsolve package, u is what we want to find, and
#p are the parameters/ numbers that we use for it
function step(state,dt)
    u_guess = guess(state,dt) #compute a guess for the solution using euler iteration
    p=[vcat(state.c,state.h),dt,state.Δx] #combine into the parameters

    #compute the sparseness of the jacobian
    jac_sparsity = ADTypes.jacobian_sparsity(
        (du, u) -> zero_func(du, u, p), similar(u_guess), u_guess, TracerSparsityDetector())
    sparse_zero_func =NonlinearFunction(zero_func; jac_prototype = jac_sparsity)
    #solve the sparse problem to make it faster
    prob = NonlinearProblem(sparse_zero_func, u_guess, p)
    sol =solve(prob,NewtonRaphson(),reltol=1e-4) #this tolerance says that we stop when Δh/h = 1e-4
    n=div(length(sol),2)
    return sol[n+1:2*n],sol[1:n] #return h and c in that order
end

#we use a step halving method to adapt the time
#compute the step once, and then with two half steps
#if the error is bigger than the res then half the time step and repeat
function two_step(state,dt)
    one = step(state,dt) #single step
    two_1 = step(state,dt/2)#double step part 1
    two = step(fornberg.state(two_1[1],two_1[2],state.coords,state.Δx,state.t+dt/2),dt/2)
    return one,two
end


function double_step_error(state,dt,res,dtmax)
    minh0 = minimum(state.h) #initial minimum thickness
    dt0=copy(dt)
    one,two = two_step(state,dt0) #compute both of them
    error = maximum(abs.(one[1]-two[1])) #error is the norm difference between the two versions
    for j in 1:5 #run this until we hit the required accuracy, don't do more than 5 loops 
        #use ricardson extrapolation to get an O(dt^2) result from the step halving method 
        new_h = 2*two[1]-one[1]
        #if error 4 times smaller than needed accept the change and increase the time step for the next step
        #also check that the thickness hasn't decreased by more than 5% in one step
        if error<res/4 && all(new_h .>0) && (minimum(two[1])/minh0 >0.95)
            dt0=min(1.1*dt0,dtmax) #do not increase the step size too much
            break #end the loop if we are way below the accuracy
        #if the error is okay then break the loop
        elseif error<res  && all(new_h .>0) && (minimum(two[1])/minh0 >0.95)
            break
        #if the error and change is too big then half the step
        elseif error > res || any(new_h .<0) || (minimum(two[1])/minh0 <0.95)
            dt0/=2
            one,two = two_step(state,dt0)
            error = maximum(abs.(one[1]-two[1]))
        end
    end
    new_h = 2*two[1]-one[1]
    new_c = 2*two[2]-one[2]
    #return the new state
    return fornberg.state(new_h,new_c,state.coords,state.Δx,state.t+dt), dt0, error
end

#performs one chunk of the calculation, i.e. iterates until the thickness 
#has changed by 5%
#data is a vector of states, dt values, and the error in the previous computation
function chunk(data,res,t_end,dxmax,L,ratio,dtmin,dtmax)
    state0=data[end][1]
    dt = data[end][2]
    min_thick_0 = minimum(state0.h) #initial thickness
    min_thick = min_thick_0
    #check that the thickness hasn't changed too much, and we haven't hit the end point, and dt isn't too small
    while (min_thick>= 0.95*min_thick_0) && (data[end][1].t <= t_end) && (data[end][2]>dtmin)
        new, dt, error = double_step_error(data[end][1],data[end][2],res,dtmax)
        push!(data,[new,dt]) #keep the new data and time step
        min_thick=minimum(new.h)
    end
    #compute the width and new dxmin
    min_coord, width, dxmin = fornberg.halfwidth(data[end][1].h,data[end][1].coords,dxmax)
    #remesh the last step if it's needed
    #i.e. mesh if dxmin<dxmax
    if dxmin < dxmax
        new_coords,new_Δx = fornberg.make_grid(min_coord,dxmin,dxmax,L,ratio,width)
        new_h = fornberg.spline(data[end][1].h,data[end][1].coords,new_coords,L)
        new_c = fornberg.spline(data[end][1].c,data[end][1].coords,new_coords,L)
        data_new=[fornberg.state(new_h,new_c,new_coords,new_Δx,data[end][1].t),dt]
        data[end]=data_new #replace the last step rather than double count it
    end
end

#This function runs many chunks till we either hit the end time t_end, or dt gets below the minimum value
function run(data,res,t_end,dxmax,L,ratio,dtmax,dtmin)
    while data[end][1].t < t_end && data[end][2]>dtmin
        chunk(data,res,t_end,dxmax,L,ratio,dtmin,dtmax)
    end
end

end


