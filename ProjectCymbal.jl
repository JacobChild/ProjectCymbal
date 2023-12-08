#=
ProjectCymbal.jl
Jacob Child
Brooke Lillie
December 5, 2023
Overview: We will be attempting to more accuratley model a cymbal (thus ProjectCymbal) by modeling the cymbal as an annular ring with Dirichlet-Robin boundary conditions. To do this a circular integral transform, Hankel transform, and Laplace transform were all taken then reversed (see pg 874 of https://www.et.byu.edu/~vps/ME505/IEM/09%2005%2018.pdf). The resultant equation models the amplitude of each r, θ, t location on the cymbal by doing a double summation over the given equations (uhat, and ubar).
Pseudocode: As this involves a double summation, it will be very calculation heavy, care will be taken to optimize speed of the code.
1. Calculate all of the eigenvalues beforehand, then make a vector of eigenfunctions for each eigenvalue (they are functions of r), so they can all be calculated quickly when called.
2. Do the inner summation uhat  for the ubar function (over indici m) multiplied by the eigenfunction over the norm.
3. Do the outer summation u for the ubar function (over indici n) 
4. Repeat for each R, θ, location on the cymbal at a given time t.
5. Repeat for each time t.
    5a. Plot the resultant function for each time t, doing both a contour plot and a surface (3D) plot.
    5b. Save in such a way so as to get an animation of the cymbal vibrating.

Notes: To save on computational time, make an integral dictionary, where the key is the integral, and the value is the result of the integral. This way, if the integral has already been calculated, it can just be looked up in the dictionary, rather than recalculated.
=#
using Roots, SpecialFunctions, QuadGK, Plots 

#functions
#find roots function 
function FindRoots(InFunc,NumofRoots)
    #finds the roots of a function
    #InFunc is the function to find the roots of
    #NumofRoots is the number of roots to find
    #returns an array of the roots
    roots = []
    i = 0
    while length(roots) != NumofRoots
        roots = find_zeros(InFunc,0.01,i)
        i+=1
    end
    return roots
end

#Integral Look up and calculate function 
function IntegralLookUp(nff, mff, IntegralDict)
    #dakey is the n and m index for the specific eigenvalue
    dakey = [nff, mff]
    # Before computing an integral
    if haskey(IntegralDict, dakey)
        # Use the cached result
        return IntegralDict[dakey]
    else
        # Compute the integral and store the result in the cache
        norm, error = quadgk(x -> eigenfunctions[nff][mff](x)^2 * x, r1, r2, rtol=1e-6)
        if error > 1e-6
            println("Warning: Integral error is large: $error")
        end
        IntegralDict[dakey] = norm
        return norm
    end
end


#Calculations 
n = 0:75
m = n
r1 = .05
r2 = 1
H = 2 #? arbitrary, no physical representation
#assuming theta0 = 0
#assuming t0 = 0
r0 = .5 #source location
S0 = .5 #source
w = 1 #? arbitrary, no physical representation
g = 1 #? arbitrary, no physical representation
res = 100
#Pseudocode pt. 1
#EigenValue Calculations from pg 526 (https://www.et.byu.edu/~vps/ME505/IEM/07%2000.pdf) Case 5, Dirichlet-Robin Boundary conditions 
function EigenValFunc(lf, nf)
    v = nf
    term1 = besselj(v, lf*r1) * ( -lf * bessely(v+1, lf*r2) + (H + v/r2) * bessely(v, lf*r2) )
    term2 = bessely(v, lf*r1) * ( -lf * besselj(v+1, lf*r2) + (H - v/r2) * besselj(v, lf*r2) )
    return term1 - term2
end
#find m # of roots (eigenvalues) at each n
eigenvalues = Vector{Vector{Float64}}(undef, length(n)) 
for (i, nf) in enumerate(n) #i is the index, nf is the value in n at that index
    EigenValFunc_fixed_n = lf -> EigenValFunc(lf, nf)
    eigenvalues[i] = FindRoots(EigenValFunc_fixed_n, length(m))
end
#EigenFunction Calculations
function EigenFunc(lf, nf)
    v = nf
    EF(rf) = besselj(v, lf*rf) / besselj(v, lf*r1) - bessely(v, lf*rf) / bessely(v, lf*r1)
    return EF
end
#find an eigenfunction for each eigenvalue
# Initialize the vector of eigenfunctions
eigenfunctions = Vector{Vector{Function}}(undef, length(n))
#Calculate the eigenfunctions for each eigenvalue
for (i, nf) in enumerate(n)
    # Initialize the vector of eigenfunctions for the current n
    eigenfunctions[i] = Vector{Function}(undef, length(m))
    for (j, lf) in enumerate(eigenvalues[i])
        # Store the eigenfunction for the current eigenvalue
        eigenfunctions[i][j] = EigenFunc(lf, nf)
    end
end

#Pseudocode pt. 2
# Initialize the integral cache
IntegralCache = Dict{Vector{Int}, Float64}() # Initialize the integral dictionary
#Functions for the pt 2 summations 
function ubar(thetaf, tf, nf, mf)
    #println("The current variables are thetaf = $thetaf, tf = $tf, nf = $nf, mf = $mf")
    return w^2 * S0 * besselj(nf,eigenvalues[nf][mf]*r0) * cos(nf*thetaf) * exp(-g*tf) * sin(sqrt(w^2 * eigenvalues[nf][mf]^2 - g^2)*tf) * H * tf / ( 2* pi * sqrt(w^2 * eigenvalues[nf][mf]^2 - g^2))
end
function an(rf, nf, mf)
    top = besselj(nf, eigenvalues[nf][mf]*rf)
    bottom = IntegralLookUp(nf, mf, IntegralCache)
    return top / bottom
end
function uhat(rf, thetaf, tf, nf; ms=m)
    #Do the summation over m of ubar and an 
    sum = 0
    for (i,mf) in enumerate(ms)
        sum += ubar(thetaf, tf, nf, i) * an(rf, nf, i)
    end
    return sum
end

#Pseudocode pt. 3
#Do the summation over n of uhat
function u(rf, thetaf, tf; ns=n)
    sum = 0
    for (i, nf) in enumerate(ns)
        if nf == 0
            sum += 1 / (2*pi) * uhat(rf, thetaf, tf, i)
        else
            sum += 1 / pi * uhat(rf, thetaf, tf, i)
        end      
    end
    return sum
end

#Pseudocode pt. 4
r = range(r1, r2, length=res)
theta = range(0, 2*pi, length=res)
t = range(0, 2, length=res)
#Functions 
function Instant(rs, thetas, tf)
    uInstantCymbal = zeros(length(rs), length(thetas))
    #intitalize x and y
    xf = zeros(length(rs), length(thetas))
    yf = zeros(length(rs), length(thetas))
    for (i, rf) in enumerate(rs)
        for (j, thetaf) in enumerate(thetas)
            uInstantCymbal[i,j] = u(rf, thetaf, tf)
            xf[i,j] = rf .* cos.(thetaf)

            yf[i,j] = rf .* sin.(thetaf)
        end
    end
    
    return uInstantCymbal
end

#Pseudocode pt. 5
# Convert the r-theta grid to the x-y space
xmesh = [rf*cos(thetaf) for rf in r, thetaf in theta]
ymesh = [rf*sin(thetaf) for rf in r, thetaf in theta]

#Output Calculations
UoverTime = [Instant(r, theta, t) for t in t]; #this is a vector of matrices, each matrix is a time step

#Plotting
t0 = surface(xmesh,ymesh,UoverTime[1], xlabel="x", ylabel="y", zlabel="u", title = "Cymbal Vibration at t = $(string(round(t[1],digits=2)))", right_margin=10*Plots.mm ,top_margin=2*Plots.mm, zlims = (-.00005, .00015), clims=(-.00005, .00015))
tstart = surface(xmesh,ymesh,UoverTime[2], xlabel="x", ylabel="y", zlabel="u", title = "Cymbal Vibration at t = $(string(round(t[2],digits=2)))", right_margin=8*Plots.mm ,top_margin=2*Plots.mm)
tstartsmall = surface(xmesh,ymesh,UoverTime[2], xlabel="x", ylabel="y", zlabel="u", title = "Cymbal Vibration at t = $(string(round(t[2],digits=2)))", right_margin=8*Plots.mm ,top_margin=2*Plots.mm, zlims=(-.0005,.0005))
tmid = surface(xmesh,ymesh,UoverTime[45], xlabel="x", ylabel="y", zlabel="u", title = "Cymbal Vibration at t = $(string(round(t[45],digits=2)))", right_margin=8*Plots.mm ,top_margin=2*Plots.mm)
tmidtop = surface(xmesh,ymesh,UoverTime[45], xlabel="x", zticks = false, title = "Cymbal Vibration at t = $(string(round(t[45],digits=2)))", right_margin=8*Plots.mm ,top_margin=2*Plots.mm,camera=(0,90))
savefig(t0, "t0.png")
savefig(tstart, "tstart.png")
savefig(tstartsmall, "tstartsmall.png")
savefig(tmid, "tmid.png")
savefig(tmidtop, "tmidtop.png")

#Create an animation using UoverTime and the top view of the cymbal
Topanim = @animate for (i, tf) in enumerate(t)
    p1 = surface(xmesh,ymesh,UoverTime[i], xlabel="x", zticks = false, title = "Cymbal Vibration at t = $(string(round(t[i],digits=2)))", right_margin=8*Plots.mm ,top_margin=2*Plots.mm,camera=(0,90))
end
gif(Topanim, "CymbalVibrationTop.gif", fps=5)



#Archive code 
#This actually runs the code 3 times for each of the plots. Iterating over UoverTime and plotting it would actually be much faster as you don't have to rerun it. I was experiencing issues with UoverTime at first, so I used this code to generate all the animations and just had it run overnight. To obtain these animations again, but much faster, just iterate over UoverTime and plot it similar to Topanim.
#= 

# Plot the results as a 3d surface over time and as a contour plot
anim = @animate for (i, tf) in enumerate(t)
    uout = [u(rf, thetaf, t[i]) for rf in r, thetaf in theta]
    p1 = surface(xmesh, ymesh, uout, xlabel="x", ylabel="y", zlabel="u", zlims=(-.0075,.0075), clims=(-.005,.005), legend=false)
    plot(p1, title = "Cymbal Vibration at t = $(string(round(tf,digits=1)))")
end
gif(anim, "CymbalVibrationHighResF.gif", fps=5)

#heatmap plot 
anim2 = @animate for (i, tf) in enumerate(t)
    uout = [u(rf, thetaf, t[i]) for rf in r, thetaf in theta]
    p2 = heatmap(uout, clims=(-.005,.005),right_margin=3*Plots.mm ,top_margin=2*Plots.mm, proj = :polar)
    plot(p2, title = "Cymbal Vibration at t = $(string(round(tf,digits=1)))")
end
gif(anim2, "CymbalVibrationHeatmapHighResF.gif", fps=5)

#contour plot
anim3 = @animate for (i, tf) in enumerate(t)
    uout = [u(rf, thetaf, t[i]) for rf in r, thetaf in theta]
    p3 = contour(r,theta,uout, xlabel="r", ylabel="theta", zlabel="u", zlims=(-.0075,.0075), clims=(-.005,.005), legend=false)
    plot(p3, title = "Cymbal Vibration at t = $(string(round(tf,digits=1)))")
end
gif(anim3, "CymbalVibrationContourHighResF.gif", fps=5)
=#