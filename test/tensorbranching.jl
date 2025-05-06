using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, truth_table, list_subtree, generate_pos_vec_local, generate_neighbors, branching_region, position_branching, ScNeighborSelector, ScRectangleSelector, tensor_branching
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


@testset "generate_pos_vec_local" begin
    n = 4
    t9_lattice = generate_TensorNQ_lattice(n)
    lattice = t9_lattice.lattice
    region = [lattice[2,j].labels[5] for j in 1:n] 

    pos0, pos1 = generate_pos_vec_local(t9_lattice, region, 1, 1)
    @test pos0 == []
    @test pos1 == [(2,1)]

    pos0, pos1 = generate_pos_vec_local(t9_lattice, region, 1, 0)
    @test pos0 == [(2,1)]
    @test pos1 == []
end


@testset "generate_neighbors for ScNeighborSelector" begin
    n = 4
    t9_lattice = generate_TensorNQ_lattice(n)
    lattice = t9_lattice.lattice
    pos_vertices =vec([lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    k = 1
    n_max = 20
    sc_target = 2
    region_selector = ScNeighborSelector(k, n_max, sc_target)

    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)
    @test length(neighbors) == length(pos_vertices)
    @test length.(neighbors) == [4,6,6,4,6,9,9,6,6,9,9,6,4,6,6,4]

    pos_vertices = pos_vertices[2:end]
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[(1,1)],Int)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)
    @test length(neighbors) == length(pos_vertices)
    @test length.(neighbors) == [5,6,4,5,8,9,6,6,9,9,6,4,6,6,4]
end


@testset "branching_region for ScNeighborSelector" begin
    n = 4
    k = 1
    n_max = 20
    sc_target = 2
    region_selector = ScNeighborSelector(k, n_max, sc_target)
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    region = branching_region(n, t9_lattice, region_selector, code, optcode)
    @test length(region) == 9
end


@testset "generate_neighbors for ScRectangleSelector" begin
    n = 4
    t9_lattice = generate_TensorNQ_lattice(n)
    lattice = t9_lattice.lattice
    pos_vertices =vec([lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    k_ud = 0
    k_lr = 2
    n_max = 20
    sc_target = 2
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target)

    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)
    @test length(neighbors) == length(pos_vertices)
    @test length.(neighbors) == [3,3,3,3,4,4,4,4,4,4,4,4,3,3,3,3]

    pos_vertices = pos_vertices[2:end]
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[(1,1)],Int)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)
    @test length(neighbors) == length(pos_vertices)
    @test length.(neighbors) == [3,3,3,3,4,4,4,3,4,4,4,3,3,3,3]
end


@testset "branching_region for ScRectangleSelector" begin
    n = 4
    k_ud = 0
    k_lr = 2
    n_max = 20
    sc_target = 2
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target)
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    region = branching_region(n, t9_lattice, region_selector, code, optcode)
    @test length(region) == 4
end


@testset "position_branching" begin
    n = 4
    solver = OptimalBranchingMIS.OptimalBranchingCore.IPSolver(optimizer=SCIP.Optimizer)
    k_ud = 0
    k_lr = 2
    n_max = 20
    sc_target = 2
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target)   
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[],Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    region_vertices = branching_region(n, t9_lattice, region_selector, code, optcode)
    branches, branch_coefficients, branch_weights = position_branching(n, t9_lattice, code, tensors, sc_target, region_vertices, solver)
    for branch in branches 
        pos1 = branch[1]
        pos0 = branch[2]
        @test length(pos1) == 1
        @test length(pos0) == 3
    end
end


@testset "tensor_branching" begin
    n = 5
    t9_lattice = generate_TensorNQ_lattice(n)
    pos1 = []
    pos0 = []
    coefficient = 1.0
    solver = OptimalBranchingMIS.OptimalBranchingCore.IPSolver(optimizer=SCIP.Optimizer)
    k_ud = 0
    k_lr = 2
    n_max = 20
    sc_target = 5
    region_selector = ScRectangleSelector(k_ud, k_lr, n_max, sc_target) 
    ccs, counting_branches = tensor_branching(n, t9_lattice, pos1, pos0, coefficient, sc_target, region_selector, solver)
    
    @test sum(counting_branches) == 10
end