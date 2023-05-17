module SymmetryFunctions

using ..Cutoff
using StaticArrays

export AbstractSymmFunction,AngularSymmFunction,RadialSymmFunction 
export RadialType2,AngularType3
export calc_one_symm_val,calc_symm_vals!,update_g_vals!

#----------------------------------------------#
#---------------Type Definitions---------------#
#----------------------------------------------#

abstract type AbstractSymmFunction{T} end 
abstract type RadialSymmFunction{T} <: AbstractSymmFunction{T} end
abstract type AngularSymmFunction{T} <: AbstractSymmFunction{T} end

#------------------------------------------------#
#-------------Symm Func Definitions--------------#
#------------------------------------------------#

struct RadialType2{T} <: RadialSymmFunction{T}
    eta::T
    r_cut::T
    type_vec::Vector
end

struct AngularType3{T} <:AngularSymmFunction{T}
    eta::T
    lambda::T
    zeta::T
    r_cut::T
    type_vec::Vector
    tpz::T
    
end

function AngularType3{T}(eta,lambda,zeta,r_cut,type_vec::Vector) where {T}
    tpz = 2.0^(1-zeta)

    return AngularType3(eta,lambda,zeta,r_cut,type_vec,tpz)
end

#------------------------------------------------------------------#
#-------------------Calculating One Symm Val-----------------------#
#------------------------------------------------------------------#
"""
    calc_one_symm_val(r2_ij,fc_ij,η)
Accepts interatomic distance squared `r2_ij`, the cutoff function 'fc_ij' and a gaussian parameter `η`  it then calculates the radial symmetry function value for a single pair of atoms.
    calc_one_symm_val(θ,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
    (position1,position2,position3,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)

returns a single symmetry function value from the double-sum. accepts `θ` the angle between ijk centred on i, and the squared distances `r2_ij`,`r2_ik`, `r2_jk`, the cutoff function values `f_ij,f_ik,f_jk` along with the symmetry funciton parameters `η`,`λ`,`ζ`, and the cutoff radius `r_cut`.

The version with `position_i` calculates the angle between positions before calculating the symmetry functions according to the previous method.  
"""
calc_one_symm_val(r2_ij,fc_ij,η) = ifelse(fc_ij!=0. && fc_ij!=1., fc_ij*exp(-η*r2_ij), 0.)

function calc_one_symm_val(θ,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
   
    g= (1+λ*θ)^ζ * exp(-η*(r2_ij+r2_ik+r2_jk)) * f_ij * f_ik * f_jk

    return g
end
function calc_one_symm_val(position1,position2,position3,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
        θ = angular_measure(position1,position2,position3,r2_ij,r2_ik)

    return calc_one_symm_val(θ,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
end

#----------------------------------------------------------------#
#------------------Total Symmetry Calculation--------------------#
#----------------------------------------------------------------#
"""
    update_g_vals!(g_vec,g_val,index1,index2,index3)
Input is the current `g_vec`tor, the `g_val`ue to update and the locations `index_i` to update. 
"""
function update_g_vals!(g_vec,g_val,index1,index2,index3)
    g_vec[index1] += g_val
    g_vec[index2] += g_val
    g_vec[index3] += g_val

    return g_vec
end
""" 
    calc_symm_vals!(positions,dist2_mat,f_mat,g_vec,symm_func::RadialType2)
Accepts `positions` for consistency with angular calculation, `dist2_mat` and `f_mat` containing the distances and cutoff functions relevant to the symmetry values, lastly accepts the symmetry function over which to iterate. `g_vec` is an N_atom vector into which the total contributions of each atom are inputted. Returns the same vector. 
"""
function calc_symm_vals!(positions,dist2_mat,f_mat,g_vec,symm_func::RadialType2)
    N=length(g_vec)
    if symm_func.type_vec == [1.,1.]
        
        for atomindex in eachindex(g_vec)
            for index2 in (atomindex+1):N
                g_val =  calc_one_symm_val(dist2_mat[atomindex,index2],f_mat[atomindex,index2],symm_func.eta)
                g_vec[atomindex] += g_val 
                g_vec[index2] += g_val
            end
        end
    else
        g_vec = zeros(N)
    end

    return g_vec
end

function calc_symm_vals!(positions,dist2_mat,f_mat,g_vec,symm_func::AngularType3)
    N = length(g_vec)
    η,λ,ζ = symm_func.eta,symm_func.lambda,symm_func.zeta
    if symm_func.type_vec == [1.,1.,1.]
        for atomindex in eachindex(g_vec)
            for index2 in (atomindex+1):N
                for index3 in (index2+1):N

                    g_val=calc_one_symm_val(positions[atomindex],positions[index2],positions[index3],dist2_mat[atomindex,index2],dist2_mat[atomindex,index3],dist2_mat[index2,index3],f_mat[atomindex,index2],f_mat[atomindex,index3],f_mat[index2,index3],η,λ,ζ)

                    g_vec = update_g_vals!(g_vec,g_val,atomindex,index2,index3)

                end
            end
        end
    else
        g_vec = zeros(N)
    end


    return symm_func.tpz*g_vec
end
#-----------------------------------------------------------------#
#--------------------Calculate Symmetry Matrix--------------------#
#-----------------------------------------------------------------#
"""
    init_symm_vecs(dist2_mat,total_symm_vec)
Prepares the symmetry matrix `g_mat` by taking the dimensions of the `dist2_mat` containing the squared distance of each atom with its pair, and `total_symm_vec` with all of the symmetry functions. 
"""
function init_symm_vecs(dist2_mat,total_symm_vec)
    g_mat=zeros(length(total_symm_vec),size(dist2_mat)[1])
    return g_mat 
end
"""
    total_symm_calc(positions,dist2_mat,f_mat,total_symm_vec)
Function to run over a vector of symmetry functions `total_symm_vec` and determining the value for each symmetry function for each atom at position `positions` with distances `dist2_mat` and a matrix of cutoff functions `f_mat` between each atom pair.
"""
function total_symm_calc(positions,dist2_mat,f_mat,total_symm_vec)
    g_mat = init_symm_vecs(dist2_mat,total_symm_vec)
    for g_index in eachindex(total_symm_vec)
        g_mat[g_index,:] = calc_symm_vals!(positions,dist2_mat,f_mat,g_mat[g_index,:],total_symm_vec[g_index])
    end
    
    return g_mat
end


end