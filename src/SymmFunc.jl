module SymmetryFunctions

using .Cutoff
using StaticArrays


abstract type AbstractSymmFunction{T} end 
abstract type RadialSymmFunction{T} <: AbstractSymmFunction{T} end
abstract type AngularSymmFunction{T} <: AbstractSymmFunction{T} end



end