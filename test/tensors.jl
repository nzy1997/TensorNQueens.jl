using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network
using OMEinsum

@testset "generate_tensor" begin
	tensor, t10, t01, t11 = generate_tensor(Float64)
	@test tensor[2, 1, 1, 1, 1, 1, 1, 1, 1] == 0.0
	@test tensor[2, 1, 1, 1, 1, 1, 1, 1, 2] == 1.0
	@test tensor[1, 1, 1, 1, 1, 1, 1, 1, 2] == 0.0
	@test tensor[1, 1, 1, 1, 1, 1, 1, 1, 1] == 1.0

	@test tensor[1, 1, 1, 1, 2, 2, 2, 2, 2] == 1.0
	@test tensor[1, 1, 1, 1, 2, 1, 2, 2, 2] == 0.0
	@test tensor[1, 2, 1, 1, 2, 2, 2, 2, 2] == 0.0
    @test tensor[1, 1, 1, 1, 1, 2, 1, 1, 1] == 0.0
	@test t10 == [1, 0]
	@test t01 == [0, 1]
	@test t11 == [1, 1]
end

@testset "TensorNQ" begin
	t = TensorNQ((1, 2, 3, 4, 5, 6, 7, 8, 9))
	@show t
end

@testset "generate_TensorNQ_lattice" begin
	tn, pos11, pos01 = generate_TensorNQ_lattice(3)
	@test tn[2, 2].labels[3] == tn[1, 1].labels[7]
	@test tn[2, 2].labels[4] == tn[1, 2].labels[6]
	@test tn[2, 2].labels[8] == tn[1, 3].labels[2]
	@test tn[2, 2].labels[1] == tn[2, 1].labels[9]
	@test tn[2, 2].labels[9] == tn[2, 3].labels[1]
	@test tn[2, 2].labels[2] == tn[3, 1].labels[8]

	@test tn[3, 1].labels[6] ∈ pos01
	@test tn[2, 2].labels[5] ∈ pos11
	@test tn[1, 3].labels[9] ∈ pos01
	@test tn[3, 3].labels[9] ∈ pos01
	@test tn[3, 3].labels[6] ∈ pos01
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