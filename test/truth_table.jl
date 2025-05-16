using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, truth_table, list_subtree, generate_pos_vec_local, generate_neighbors, branching_region, position_branching, ScNeighborSelector, ScRectangleSelector, ScRectangleD4Selector,tensor_branching, naive_tensor_branching, random_naive_tensor_branching, sc_score_weights, ExactScScorer, TruncateBondScorer
using OMEinsum
using OptimalBranching
using SCIP
using Graphs


@testset "generate_truth_table" begin
    n = 4
    t9_lattice = generate_TensorNQ_lattice(n)
    lattice = t9_lattice.lattice
    
    region_rows = [[lattice[i,j].labels[5] for j in 1:n] for i in 1:n]
    region_cols = [[lattice[i,j].labels[5] for j in 1:n] for i in 1:n]
    region_diagonals = [[lattice[i,i].labels[5] for i in 1:n], [lattice[i,n-i+1].labels[5] for i in 1:n]]
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[], Int)

    for region in region_rows
        configs = truth_table(code, tensors, region)
        @test length(configs) == n
    end
    for region in region_cols
        configs = truth_table(code, tensors, region)
        @test length(configs) == n
    end
    for region in region_diagonals
        configs = truth_table(code, tensors, region)
        @test length(configs) == n+1
    end
end
