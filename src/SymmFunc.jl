module SymmetryFunctions

using ..Cutoff
using StaticArrays

export AbstractSymmFunction,AngularSymmFunction,RadialSymmFunction 
export RadialType2,AngularType3
export calc_symmetry_function

#----------------------------------------------#
#---------------Type Definitions---------------#
#----------------------------------------------#

abstract type AbstractSymmFunction{T} end 
abstract type RadialSymmFunction{T} <: AbstractSymmFunction{T} end
abstract type AngularSymmFunction{T} <: AbstractSymmFunction{T} end

#------------------------------------------------#
#--------------Type 2 symm function--------------#
#------------------------------------------------#

struct RadialType2{T} <: RadialSymmFunction{T}
    eta::T
    r_cut::T
    type_vec::Vector
end
"""
    calc_one_symm_function(r2_ij,eta,r_cut)
Accepts interatomic distance squared `r2_ij`, gaussian parameter `η` and cutoff distance `r_cut` it then calculates the radial symmetry function value for a single pair of atoms.
"""
function calc_one_symm_function(r2_ij,η,r_cut)
    if r2_ij != 0 && r2_ij < r_cut^2 
        g_ij = exp(-η*(r2_ij)) * cutoff_function( sqrt(r2_ij)/r_cut )
        return g_ij
    else
        return 0.
    end
end
"""
    calc_symm_function(dist2_matrix,index,symmfunc::RadialType2)
calculate the total symmetry function value `G` centred on atom `index` by iterating over the interatomic distances in `dis2_mat` using the symmetry parameters contained in `symmfunc`. 
"""
function calc_symm_function(dist2_matrix,index,symmfunc::RadialType2)
    eta = symmfunc.eta
    r_cut = symmfunc.r_cut    
    gvec = calc_one_symm_function.(dist2_matrix[:,index],eta,r_cut)
    return sum(gvec)   
end


#------------------------------------------------#
#--------------Type 3 symm function--------------#
#------------------------------------------------#
struct AngularType3{T} <:AngularSymmFunction{T}
    eta::T
    lambda::T
    zeta::T
    r_cut::T
    type_vec::Vector
end
"""
    calc_one_symm_value(θ,r2_ij,r2_ik,r2_jk,r_cut,η,λ,ζ)
returns a single symmetry function value from the double-sum. accepts `θ` the angle between ijk centred on i, and the squared distances `r2_ij`,`r2_ik`, `r2_jk` along with the symmetry funciton parameters `η`,`λ`,`ζ`, and the cutoff radius `r_cut`.
"""
function calc_one_symm_value(θ,r2_ij,r2_ik,r2_jk,r_cut,η,λ,ζ)
    d_cut=r_cut^2
    if r2_ij > d_cut || r2_ik > d_cut || r2_jk > d_cut 
        g = 0.
    else
        g = (1+λ*θ)^ζ*exp(-η*(r2_ij+r2_ik+r2_jk))*cutoff_function(sqrt(r2_ij)/r_cut)*cutoff_function(sqrt(r2_ik)/r_cut)*cutoff_function(sqrt(r2_jk)/r_cut)
    end
    
    return g
end
"""
    calc_one_symm_func(position1,position2,position3,symmfunc)
    calc_one_symm_func(position1,position2,position3,r2_ij,r2_ik,r2_jk,symmfunc)
function to calculate the symmetry function value with parameters in `symmfunc` of a single triplet of atoms at `position1` `position2` `position3` centred on position 1. The first method calculates the interatomic distances `r2_ij`,`r2_ik`,`r2_jk` and cosθ labelled as `θ`, returning `g` and the distances, the second accepts these as arguments and only returns `g`.
"""
function calc_one_symm_func(position1,position2,position3,symmfunc::AngularType3)
    θ,r2_ij,r2_ik,r2_jk = angular_measure(position1,position2,position3)
    g = calc_one_symm_value(θ,r2_ij,r2_ik,r2_jk,symmfunc.r_cut,symmfunc.eta,symmfunc.lambda,symmfunc.zeta)
    return g,r2_ij,r2_ik,r2_jk
end
function calc_one_symm_func(position1,position2,position3,r2_ij,r2_ik,r2_jk,symmfunc::AngularType3)
    θ = angular_measure(position1,position2,position3,r2_ij,r2_ik)
    g = calc_one_symm_value(θ,r2_ij,r2_ik,r2_jk,symmfunc.r_cut,symmfunc.eta,symmfunc.lambda,symmfunc.zeta)
    return g
end
"""
    calc_symmetry_function(positions,dis2_mat,index,symmfunc::AngularType3)
accepts a vector of `positions` corresponding to a configuration,   `dis2_mat` containing the interatomic distances, `index` corresponding to the central atom and `symmfunc` with the parameters in use. 
"""
function calc_symm_function(positions,dis2_mat,index,symmfunc::AngularType3)
    g_vec = [calc_one_symm_func(positions[index],positions[j],positions[k],dis2_mat[index,j],dis2_mat[index,k],dis2_mat[j,k],symmfunc) for j=(1:55) if j !=index for k= 1:j-1 if k != index]
    g = 2^(1-symmfunc.zeta)*sum(g_vec)
    return g
end

#--------------------------------------------------------#
#-------------------Combined function--------------------#
#--------------------------------------------------------#\
"""
    calc_symmetry_function(positions,dis2_mat,index,symmfunc::RadialType2)
    calc_symmetry_function(positions,dis2_mat,index,symmfunc::AngularType3)
    calc_symmetry_function(positions,dis2_mat,index,symmfunc,g_min,g_max)
End-use function for calculating symmetry functions, designed as a curry-function to ensure the correct method is called with unchanged input. accepts a vector `positions` a matrix `dis2_mat` of square distances, the `index` of the atom about which we calculate and the `symmfunc` containing the parameters. 

    -optional inclusion of `g_max` `g_min` designed to rescale the symmetry function to ensure g∈[0,1] 

"""
function calc_symmetry_function(positions,dis2_mat,index,symmfunc::RadialType2)
    if symmfunc.type_vec == [1.,1.]
        g = calc_symm_function(dis2_mat,index,symmfunc)
    else
        g=0.
    end

    return g
end
function calc_symmetry_function(positions,dis2_mat,index,symmfunc::AngularType3)
    if symmfunc.type_vec == [1.,1.,1.]
        g = calc_symm_function(positions,dis2_mat,index,symmfunc::AngularType3)
    else 
        g = 0.
    end

    return g
end
function calc_symmetry_function(positions,dis2_mat,index,symmfunc,g_min,g_max)
    g_unscaled = calc_symmetry_function(positions,dis2_mat,index,symmfunc)
    return (g_unscaled - g_min)/(g_max - g_min) 
end
function calc_symmetry_function(positions,dis2_mat,index,symmfunc,scaledata::Vector)
    return calc_symmetry_function(positions,dis2_mat,index,symmfunc,scaledata[1],scaledata[2])
end

end