module Cutoff

using LinearAlgebra
using StaticArrays

export distance2,get_distance2_mat
export cutoff_function


distance2(a,b) = (a-b)⋅(a-b)

get_distance2_mat(pos) = [distance2(a,b) for a in pos, b in pos]

function cutoff_function(r_scaled)
    
    cutoff= 0.5*(cos(π*r_scaled) + 1)
    

    return cutoff
end
cutoff_function(r_ij,r_cut) = cutoff_function(r_ij/r_cut)


end