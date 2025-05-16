using Test
using TensorNQueens
using TensorNQueens: set_logging, generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, truth_table, list_subtree, generate_pos_vec_local, generate_neighbors, branching_region, position_branching, ScNeighborSelector, ScRectangleSelector, ScRectangleD4Selector,tensor_branching, naive_tensor_branching, random_naive_tensor_branching, sc_score_weights, ExactScScorer, TruncateBondScorer
using OMEinsum
using OptimalBranching
using SCIP
using Graphs


@testset "position_branching for TruncateBondScorer" begin
    n = 4
    solver = OptimalBranchingMIS.OptimalBranchingCore.IPSolver(optimizer=SCIP.Optimizer)
    k_ud = 0
    k_lr = 2
    n_max = 20
    sc_target = 2
    bond_limit = 4
    scorer = TruncateBondScorer(bond_limit,true,false)
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target)   
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    region_vertices = branching_region(n, t9_lattice, region_selector, code, optcode)
    branches, branch_coefficients, branch_weights = position_branching(n, t9_lattice, code, tensors, sc_target, region_vertices, solver, scorer)
    @test length(branches) == 2
    @test Set([length(branches[1][1]),length(branches[2][1])]) == Set([0,1])
    @test Set([length(branches[1][2]),length(branches[2][2])]) == Set([1,3])
end


@testset "tensor_branching for SC_based region selector and D4_based scorer" begin
    n = 10
    t9_lattice = generate_TensorNQ_lattice(n)
    pos1 = []
    pos0 = []
    coefficient = 1.0
    solver = OptimalBranchingMIS.OptimalBranchingCore.IPSolver(optimizer=SCIP.Optimizer)
    k_ud = 0
    k_lr = 5
    n_max = 20
    sc_target = 20
    bond_limit = 4
    scorer = TruncateBondScorer(bond_limit, true, false)
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target) 
    output_file = "test_log/n=$(n)_sc_target=$(sc_target)_SCregion_D4scorer_rec=$(k_lr)_$(k_ud).log"
    set_logging(true, output_file)
    ccs, counting_branches = tensor_branching(n, t9_lattice, pos1, pos0, coefficient, sc_target, region_selector, solver, scorer, 1)
    
    @test abs(sum(counting_branches) - 724) < 1e-6
end


@testset "tensor_branching for D4_based region selector and D4_based scorer" begin
    n = 10
    t9_lattice = generate_TensorNQ_lattice(n)
    pos1 = []
    pos0 = []
    coefficient = 1.0
    solver = OptimalBranchingMIS.OptimalBranchingCore.IPSolver(optimizer=SCIP.Optimizer)
    k_ud = 0
    k_lr = 5
    n_max = 20
    sc_target = 20
    bond_limit = 4
    scorer = TruncateBondScorer(bond_limit,true,false)
    region_selector = ScRectangleD4Selector(k_ud, k_lr, n_max, sc_target) 
    output_file = "test_log/n=$(n)_sc_target=$(sc_target)_D4region_D4scorer_rec=$(k_lr)_$(k_ud).log"
    set_logging(true, output_file)
    distances, scs, counting_branches = tensor_branching(n, t9_lattice, pos1, pos0, coefficient, sc_target, region_selector, solver, scorer, 1)
    
    @test abs(sum(counting_branches) - 724) < 1e-6
end


@testset "naive_tensor_branching" begin
    n = 10
    t9_lattice = generate_TensorNQ_lattice(n)
    pos1 = []
    pos0 = []
    sc_target = 20
    output_file = "test_log/n=$(n)_sc_target=$(sc_target)_naive_row.log"
    set_logging(true, output_file)
    ccs, counting_branches = naive_tensor_branching(n, t9_lattice, pos1, pos0, sc_target, 1)
    
    @test abs(sum(counting_branches) - 724) < 1e-6
end
