# MCMachineLearningPotentialInterface

# About the Project
A package developed in conjunction with the RuNNer2.0 Julia Interface (https://gitlab.com/runner-suite/runner-julia-interface.git) to calculate the symmetry functions and energies of a neural network that has been trained using the RuNNer program. 
This was designed as a complete interface for use with the ParallelTemperingMonteCarlo package (https://github.com/ElkePahl/ParallelTemperingMonteCarlo.jl.git) in order to run atomic simulations using MLP's. 

# Usage
- Requires a compiled librunnerjulia.so file from the RuNNer2.0 Julia Interface
- Requires the symmfunctions.out, scaling.data and weights.XXX.data file from a RuNNer-trained NNP
- examples of use are in the ParallelTemperingMonteCarlo package (TODO: add basic script showing use)

# Contact
Gray Hunter: gray.hunter@auckland.ac.nz





[![Build Status](https://github.com/grayseff/MachineLearningPotential.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/grayseff/MachineLearningPotential.jl/actions/workflows/CI.yml?query=branch%3Amain)
