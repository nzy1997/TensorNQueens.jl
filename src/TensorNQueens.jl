module TensorNQueens

using OMEinsum
using BitBasis
using TropicalSweepContractor

using TropicalSweepContractor: PlanarTensor

include("tensors.jl")
include("tensor8.jl")
include("tensor3.jl")
include("sweepcontract.jl")
end
