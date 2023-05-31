module SymmetryFunctions

using ..Cutoff
using StaticArrays

export AbstractSymmFunction,AngularSymmFunction,RadialSymmFunction 
export RadialType2,AngularType3
export calc_one_symm_val,calc_symm_vals!,update_g_vals!
export total_symm_calc

#------------------------------------------------#
#----------------Type Definitions----------------#
#------------------------------------------------#

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
    G_offset::T
    G_norm::T
end
function RadialType2{T}(eta,r_cut,type_vector) where {T}
    return RadialType2(eta,r_cut,type_vector,0.,1.)
end
function RadialType2{T}(eta,r_cut,type_vector::Vector,G_vals::Vector) where {T}
    G_norm = 1/(G_vals[1] - G_vals[2])
    G_offset = -G_vals[1]*G_norm
    return RadialType2(eta,r_cut,type_vector,G_offset,G_norm)
end

struct AngularType3{T} <:AngularSymmFunction{T}
    eta::T
    lambda::T
    zeta::T
    r_cut::T
    type_vec::Vector
    tpz::T
    G_offset::T
    #G_norm::T
    
end

function AngularType3{T}(eta,lambda,zeta,r_cut,type_vec::Vector) where {T}
    tpz = 2.0^(1-zeta)
    return AngularType3(eta,lambda,zeta,r_cut,type_vec,tpz,0.)
end
function AngularType3{T}(eta,lambda,zeta,r_cut,type_vector::Vector,G_vals::Vector) where {T}
    G_norm = 1/(G_vals[1] - G_vals[2])
    G_offset = -G_vals[1]*G_norm
    tpz = 2.0^(1-zeta)*G_offset
    return AngularType3(eta,lambda,zeta,r_cut,type_vector,tpz,G_offset)
end
#------------------------------------------------------------------#
#-------------------Calculating One Symm Val-----------------------#
#------------------------------------------------------------------#
"""
    exponential_part(η,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk)
calculates the exponential portion of the symmetry function for the angular symmetry function. Preserves the values we can maintain throughout iterating over theta
"""
exponential_part(η,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk) = exp(-η*(r2_ij+r2_ik+r2_jk))* f_ij * f_ik * f_jk
"""
    theta_part(θ,λ,ζ)
Calculates the angular portion of a single symmetry function, this requires iteration over each of the three angles.
"""
theta_part(θ,λ,ζ) = (1+λ*θ)^ζ
"""
    symmfunc_calc(θ_vec,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
Calculates the three g_values corresponding to the three atoms iterated over, builds the foundation of the total symm function as calculated below.
"""
function symmfunc_calc(θ_vec,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)

    exp_part = exponential_part(η,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk)
    g_values = [exp_part* theta_part(θ,λ,ζ) for θ in θ_vec]
    
    return g_values
end
"""
    update_g_vals!(g_vec,g_vals,atomindex,index2,index3)
function to correctly update the symmvalues 'g_vals' at the indices in 'g_vec'
 """
function update_g_vals!(g_vec,g_vals,atomindex,index2,index3)

    g_vec[atomindex] += g_vals[1]
    g_vec[index2] += g_vals[2]
    g_vec[index3] += g_vals[3]
    
    return g_vec
end

"""
    calc_one_symm_val(r2_ij,fc_ij,η)
Accepts interatomic distance squared `r2_ij`, the cutoff function 'fc_ij' and a gaussian parameter `η`  it then calculates the radial symmetry function value for a single pair of atoms.
    calc_one_symm_val(θ,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
    (position1,position2,position3,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)

returns a single symmetry function value from the double-sum. accepts `θ` the angle between ijk centred on i, and the squared distances `r2_ij`,`r2_ik`, `r2_jk`, the cutoff function values `f_ij,f_ik,f_jk` along with the symmetry funciton parameters `η`,`λ`,`ζ`, and the cutoff radius `r_cut`.

The version with `position_i` calculates the angle between positions before calculating the symmetry functions according to the previous method.  
"""
calc_one_symm_val(r2_ij,fc_ij,η) = ifelse(fc_ij!=0. && fc_ij!=1., fc_ij*exp(-η*r2_ij), 0.)


function calc_one_symm_val(position1,position2,position3,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)
    θ_vec = all_angular_measure(position1,position2,position3,r2_ij,r2_ik,r2_jk)
    
    g_vals = symmfunc_calc(θ_vec,r2_ij,r2_ik,r2_jk,f_ij,f_ik,f_jk,η,λ,ζ)

    return g_vals
end

#----------------------------------------------------------------#
#------------------Total Symmetry Calculation--------------------#
#----------------------------------------------------------------#

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
    if symm_func.type_vec == [1.,1.,1.]
        η,λ,ζ = symm_func.eta,symm_func.lambda,symm_func.zeta
        for atomindex in eachindex(g_vec)
            for index2 in (atomindex+1):N
                for index3 in (index2+1):N

                    g_vals=calc_one_symm_val(positions[atomindex],positions[index2],positions[index3],dist2_mat[atomindex,index2],dist2_mat[atomindex,index3],dist2_mat[index2,index3],f_mat[atomindex,index2],f_mat[atomindex,index3],f_mat[index2,index3],η,λ,ζ)

                    g_vals .*= symm_func.tpz 
                    
                    g_vec = update_g_vals!(g_vec,g_vals,atomindex,index2,index3)
                end
            end
        end
        
    else
        g_vec = zeros(N)
    end


    return g_vec

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