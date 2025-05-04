using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network,generate_masked_8_tensor_network
using OMEinsum

@testset "generate_tensor" begin
	tensor = generate_tensor(Float64)
	@test tensor[2, 1, 1, 1, 1, 1, 1, 1, 1] == 0.0
	@test tensor[2, 1, 1, 1, 1, 1, 1, 1, 2] == 1.0
	@test tensor[1, 1, 1, 1, 1, 1, 1, 1, 2] == 0.0
	@test tensor[1, 1, 1, 1, 1, 1, 1, 1, 1] == 1.0

	@test tensor[1, 1, 1, 1, 2, 2, 2, 2, 2] == 1.0
	@test tensor[1, 1, 1, 1, 2, 1, 2, 2, 2] == 0.0
	@test tensor[1, 2, 1, 1, 2, 2, 2, 2, 2] == 0.0
    @test tensor[1, 1, 1, 1, 1, 2, 1, 1, 1] == 0.0
end

@testset "generate_TensorNQ_lattice" begin
	tn_lattice = generate_TensorNQ_lattice(3)
	@test tn_lattice.lattice[2, 2].labels[3] == tn_lattice.lattice[1, 1].labels[7]
	@test tn_lattice.lattice[2, 2].labels[4] == tn_lattice.lattice[1, 2].labels[6]
	@test tn_lattice.lattice[2, 2].labels[8] == tn_lattice.lattice[1, 3].labels[2]
	@test tn_lattice.lattice[2, 2].labels[1] == tn_lattice.lattice[2, 1].labels[9]
	@test tn_lattice.lattice[2, 2].labels[9] == tn_lattice.lattice[2, 3].labels[1]
	@test tn_lattice.lattice[2, 2].labels[2] == tn_lattice.lattice[3, 1].labels[8]

	@test tn_lattice.lattice[3, 1].labels[6] ∈ tn_lattice.pos01
	@test tn_lattice.lattice[2, 2].labels[5] ∈ tn_lattice.pos11
	@test tn_lattice.lattice[1, 3].labels[9] ∈ tn_lattice.pos01
	@test tn_lattice.lattice[3, 3].labels[9] ∈ tn_lattice.pos01
	@test tn_lattice.lattice[3, 3].labels[6] ∈ tn_lattice.pos01
end

@testset "generate_tensor_network" begin
    n = 7
    code, tensors = generate_tensor_network(n, Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    @test optcode(tensors...)[] == 40
end

@testset "generate_8_tensor_network" begin
    n = 7
    code, tensors = generate_8_tensor_network(n, Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    @test optcode(tensors...)[] == 40
end

@testset "generate_3_tensor" begin
    t = TensorNQueens.generate_3_tensor(Int)
    t8 = TensorNQueens.generate_8_tensor(Int)
    @test ein"abi,cdi,efi,ghi->aceghfdb"(t,t,t,t) == t8
end

@testset "generate_3_tensor_network" begin
    n = 7
    code, tensors = generate_3_tensor_network(n, Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    @info contraction_complexity(optcode, uniformsize(optcode, 2))
    @test optcode(tensors...)[] == 40
end

@testset "benchmarking" begin
    # for n in 1:15
    for n in 1:1
        code, tensors = generate_3_tensor_network(n, Int)
        time_start = time()
        @info "n = $n"
        optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
        @info contraction_complexity(optcode, uniformsize(optcode, 2))
        time_end1 = time()
        @info optcode(tensors...)[]
        time_end2 = time()
        @info "time = $(time_end1 - time_start), $(time_end2 - time_end1)"
        println()
    end
end

@testset "generate_pos_vec" begin
    n = 4
    pos0, pos1 = TensorNQueens.generate_pos_vec(n, 0b011010111, 0b001000111)
    @test pos0 == [(1, 2), (4, 2)]
    @test pos1 == [(1, 1), (2, 1), (3, 1), (3, 2)]
end

@testset "generate_masked_3_tensor_network" begin
    n = 5
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[(2,1)],[], Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    ans1 = optcode(tensors...)[]

    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[],[(2,1)], Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    ans2 = optcode(tensors...)[]
    @test ans1 + ans2 == 10
end


@testset "generate_masked_3_tensor_network" begin
    n = 5
    t9_lattice = generate_TensorNQ_lattice(n)
    mask = 0b11
    val_vec = [0b00,0b01,0b10]
    s = 0 
    for val in val_vec
        code, tensors = generate_masked_3_tensor_network(n,t9_lattice,mask,val, Int)
        optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
        s += optcode(tensors...)[]
        @show s
    end
    @test s == 10
end

@testset "benchmarking on branching" begin
    n = 9
    col = 5
    s = Int[]
    for row in 1:5
        t9_lattice = generate_TensorNQ_lattice(n)
        code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[(row,col)],[], Int)
        optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
        @info "row = $row"
        @info contraction_complexity(optcode, uniformsize(optcode, 2))
        s1 = optcode(tensors...)[]
        @info "s1 = $s1"
        push!(s,s1)
    end
    @test 2 * sum(s[1:4])+s[5] == 352
end


@testset "generate_masked_3_tensor_network" begin
    n = 5
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,[(1,3)],[], Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    ans1 = optcode(tensors...)[]
    @test ans1 == 2
end

@testset "generate_masked_8_tensor_network" begin
    solver = TreeSA()
    # for n in 5:5
    for n in 5:28
        code, tensors = generate_8_tensor_network(n, Int)
        t9_lattice = generate_TensorNQ_lattice(n)

        code2, tensors2 = generate_masked_8_tensor_network(n,t9_lattice,[(1,n÷2 +1)],[], Int)
        time_start = time()
        @info "n = $n"
        optcode = optimize_code(code, uniformsize(code, 2), solver)
        @info contraction_complexity(optcode, uniformsize(optcode, 2))
        time_end1 = time()

        optcode2 = optimize_code(code2, uniformsize(code2, 2), solver)
        @info contraction_complexity(optcode2, uniformsize(optcode2, 2))
        time_end2 = time()
        @info "time = $(time_end1 - time_start), $(time_end2 - time_end1)"
        println()
    end
end

@testset "generate_masked_8_tensor_network" begin
    n = 5
    t9_lattice = generate_TensorNQ_lattice(n)
    code, tensors = generate_masked_8_tensor_network(n,t9_lattice,[(1,3)],[], Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    ans1 = optcode(tensors...)[]
    @test ans1 == 2
end

@testset "generate_masked_8_tensor_network" begin
    n = 5
    t9_lattice = generate_TensorNQ_lattice(n)
    mask = 0b11
    val_vec = [0b00,0b01,0b10]
    s = 0 
    for val in val_vec
        code, tensors = generate_masked_8_tensor_network(n,t9_lattice,mask,val, Int)
        optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
        s += optcode(tensors...)[]
        @show s
    end
    @test s == 10
end