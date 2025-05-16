using TensorNQueens
using Test

@testset "TensorNQueens.jl" begin
    include("tensors.jl")
end

@testset "TensorNQueens.jl" begin
    include("tensorbranching.jl")
end

@testset "TensorNQueens.jl" begin
    include("truth_table.jl")
end

@testset "TensorNQueens.jl" begin
    include("branching_region.jl")
end