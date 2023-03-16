module Cutoff

using LinearAlgebra
using StaticArrays

export distance2,get_distance2_mat
export cutoff_function

#----------------------------------------------------------------#
#-------------------------Measurements---------------------------#
#----------------------------------------------------------------#

#Beginning with basic distance funtions required throughout MLP calculations
"""
    distance2(a,b)
squared distance of two vectors `a` `b` 
"""
distance2(a,b) = (a-b)⋅(a-b)
"""
    get_distance2_mat(pos)
given a vector called `pos` comprised of (ideally) static vectors we return a lengthXlength symmetric matrix of the squared distance
"""
get_distance2_mat(pos) = [distance2(a,b) for a in pos, b in pos]

#------------------------------------------------------------------#
#----------------------Type 2 cutoff function----------------------#
#------------------------------------------------------------------#
"""
    cutoff_function(r_scaled)
    cutoff_function(r_ij,r_cut)
Implementation of the type 2 cutoff function. Either accepts scaled radius `r_scaled` or the interatomiic distance `r_ij` and the cutoff radius `r_cut`. Calculation is described in the RuNNer documentation, given as 1/2 (cos(πx) + 1) where x is (r_ij - r_i,c)/(rc - r_i,c). As an inner cutoff is not used by the potentials we are interested in, we have not included a method. 
"""
function cutoff_function(r_scaled)
    
    cutoff= 0.5*(cos(π*r_scaled) + 1)
    

    return cutoff
end
cutoff_function(r_ij,r_cut) = cutoff_function(r_ij/r_cut)


end