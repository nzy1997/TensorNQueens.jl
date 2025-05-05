using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, truth_table, list_subtree, generate_pos_vec_local, generate_neighbors, branching_region, position_branching, ScNeighborSelector, ScRectangleSelector, tensor_branching
using OMEinsum
using OptimalBranching
using SCIP
using Graphs

# @testset "generate_masked_3_tensor_network" begin
#     niters = 1000
#     for n in 28:28
#         code, tensors = generate_3_tensor_network(n, Int)
#         t9_lattice = generate_TensorNQ_lattice(n)

#         code2, tensors2 = generate_masked_3_tensor_network(n,t9_lattice,[(1,n÷2 +1)],[], Int)
        
#         time_start = time()
#         @info "n = $n"
#         optcode = optimize_code(code, uniformsize(code, 2), TreeSA(niters = niters))
#         @info contraction_complexity(optcode, uniformsize(optcode, 2))
#         time_end1 = time()

#         optcode2 = optimize_code(code2, uniformsize(code2, 2), TreeSA(niters = niters))
#         @info contraction_complexity(optcode2, uniformsize(optcode2, 2))
#         time_end2 = time()
#         @info "time = $(time_end1 - time_start), $(time_end2 - time_end1)"
#         println()
#     end
# end

@testset "tensor_branching" begin
    n = 12
    t9_lattice = generate_TensorNQ_lattice(n)
    # code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    # optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    # cc = contraction_complexity(optcode, uniformsize(optcode, 2))
    # println(cc)
    pos1 = []
    pos0 = []
    coefficient = 1.0
    solver = OptimalBranchingMIS.OptimalBranchingCore.IPSolver(optimizer=SCIP.Optimizer)
    k_ud = 0
    k_lr = 3
    n_max = 20
    sc_target = 25
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target) 
    ccs, counting_branches = tensor_branching(n, t9_lattice, pos1, pos0, coefficient, sc_target, region_selector, solver)
    println(ccs)
    println(counting_branches)
    println(sum(counting_branches))
end