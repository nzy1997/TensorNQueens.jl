using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_8_tensor_network
using OMEinsum

@testset "generate_masked_3_tensor_network" begin
    solver = TreeSA(niters = 500)
    for n in 28:28
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