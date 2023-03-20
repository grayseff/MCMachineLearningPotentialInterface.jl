module SymmetryFunctions

using ..Cutoff
using StaticArrays

export AbstractSymmFunction,AngularSymmFunction,RadialSymmFunction 
export RadialType2,AngularType3

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
end


struct AngularType3{T} <:AngularSymmFunction{T}
    eta::T
    r_cut::T
    lambda::T
    zeta::T
end


end