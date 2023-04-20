module Cutoff

using LinearAlgebra
using StaticArrays

export distance2,get_distance2_mat,angular_measure
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
"""
    angular_measure(a,b,c,r2ij,r2ik)
    angular_measure(a,b,c)
angular measure accepts three vectors `a`,`b`,`c` and can either accept or calculate the squared distances between them `r2_ab`,`r2_bc`, centred on vector `a`. Returns cos(θ) labelled as θ: the angular measure.
"""
function angular_measure(a,b,c,r2ab,r2ac)
    θ = (a - b)⋅(a - c)/sqrt(r2ab*r2ac) 
    return θ
end
function angular_measure(a,b,c)
    r2_ab,r2_ac,r2_bc = distance2(a,b),distance2(a,c),distance2(b,c)
    θ = angular_measure(a,b,c,r2_ab,r2_ac)   
    return θ,r2_ab,r2_ac,r2_bc
end
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
cutoff_function(r_ij,r_cut) = ifelse(r_ij<r_cut ,cutoff_function(r_ij/r_cut),0.)


end