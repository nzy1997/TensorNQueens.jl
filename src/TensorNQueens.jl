module TensorNQueens

using OMEinsum
using BitBasis
using OptimalBranching
using OptimalBranching.OptimalBranchingCore, OptimalBranching.OptimalBranchingMIS
using AbstractTrees, TreeWidthSolver
using JLD2
using Graphs
using SCIP
using Graphs

export ScNeighborSelector

include("utils.jl")
include("tensors.jl")
include("tensor8.jl")
include("tensor3.jl")
include("tensorbranching.jl")
end
