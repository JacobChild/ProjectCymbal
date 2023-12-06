#=
ProjectCymbal.jl
Jacob Child
Brooke Lillie
December 5, 2023
Pseudocode: Calculate all of the eigenvalues beforehand, then make a vector of eigenfunctions for each eigenvalue (they are functions of r). 
=#
using Roots, SpecialFunctions 

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


#Calculations 
n = [0,1,2]
m = [0,1,2]
r1 = .05
r2 = 1
H = 2
#assuming theta0 = 0
#assuming t0 = 0
r0 = .5
S0 = 1 #source
w = 1 #?

#EigenValue Calculations 
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

for (i, nf) in enumerate(n)
    # Initialize the vector of eigenfunctions for the current n
    eigenfunctions[i] = Vector{Function}(undef, length(m))
    for (j, lf) in enumerate(eigenvalues[i])
        # Store the eigenfunction for the current eigenvalue
        eigenfunctions[i][j] = EigenFunc(lf, nf)
    end
end

#Functions for the summations 
function ubar(thetaf, tf, nf, mf)
    return w^2 * S0 * besselj(nf,eigenvalues[nf][mf]*r0) * cos(nf*thetaf) * exp(-g*tf) * sin(sqrt(w^2 * eigenvalues[nf][mf]^2 - g^2)*tf) * H * tf / ( 2* pi * sqrt(w^2 * eigenvalues[nf][mf]^2 - g^2))
end