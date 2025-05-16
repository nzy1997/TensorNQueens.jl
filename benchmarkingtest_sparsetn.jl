using OMSparseEinsum
using Test
using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network,generate_masked_8_tensor_network, find_max_in_file
using OMEinsum
using Dates
using OMEinsumContractionOrders
using JuMP
using HiGHS
using LinearAlgebra


function find_max_in_file(filename)
    max_value = -Inf
    open(filename, "r") do file
        for line in eachline(file)
            current_value = parse(Float64, line)
            if current_value > max_value
                max_value = current_value
            end
        end
    end
    
    return max_value
end

function solve_nqueens_IP(N)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:N, 1:N], Bin)
    for i in 1:N
        @constraint(model, sum(x[i, :]) == 1)
        @constraint(model, sum(x[:, i]) == 1)
    end
    for i in -(N - 1):(N-1)
        @constraint(model, sum(LinearAlgebra.diag(x, i)) <= 1)
        @constraint(model, sum(LinearAlgebra.diag(reverse(x; dims = 1), i)) <= 1)
    end
    optimize!(model)
    assert_is_solved_and_feasible(model)
    solution = round.(Int, value.(x))
   
    return solution
end


function fix_one_kind_encirclement(n::Int,L::Int)
    solution = solve_nqueens_IP(2*L)
    pos0 = []
    pos1 = []
    for i in 1:2*L
        for j in 1:2*L
            i_pos = i
            j_pos = j
            if i_pos > L
                i_pos = i_pos - 2*L + n
            end
            if j_pos > L
                j_pos = j_pos - 2*L + n
            end
            if solution[i,j] == 1
                push!(pos1, (i_pos,j_pos))
            else
                push!(pos0, (i_pos,j_pos))
            end
        end
    end
    return pos1, pos0
end


@testset "test tc and sc for sparse tensor network" begin
    for n in 4:13
        code, dense = generate_tensor_network(n, Int)
        t1 = @elapsed tensors = map(t -> SparseTensor(t), dense)
        t2 = @elapsed size_dict = OMEinsum.get_size_dict(code.ixs, tensors)
        t3 = @elapsed optcode = optimize_code(code, size_dict, TreeSA())

        contraction_stats = @timed value = optcode(tensors...)
        contraction_time_taken = contraction_stats.time  # seconds
        contraction_memory_bytes = contraction_stats.bytes  # bytes
        contraction_memory_mb = contraction_memory_bytes / (1024 * 1024)
        sc = log2(contraction_memory_mb) + 17.0

        cc = OMEinsumContractionOrders.contraction_complexity(optcode, size_dict)

        open("data/profile/sparse_tn.txt", "a") do io
            println(io, "n:$n", ", value:",value[1])
            println(io, "Time for sparse map: $(round(t1, digits=6)) seconds")
            println(io, "Time for size_dict: $(round(t2, digits=6)) seconds")
            println(io, "Time for optimize_code: $(round(t3, digits=6)) seconds")
            println(io, "Time for contraction: $(round(contraction_time_taken, digits=6)) seconds")
            println(io, "Memory for contraction: $(round(contraction_memory_mb, digits=6)) MB")
            println(io, "SC: $(round(sc, digits=6))")
            println(io, "contraction complexity: ")
            println(io, cc)
            println(io, "--------------------------------")
        end
    end
end


@testset "test num of non-zero elements in the contraction process for sparse tensor network" begin
    for n in 4:13
        code, dense = generate_tensor_network(n, Int)
        tensors = map(t -> SparseTensor(t), dense)
        size_dict = OMEinsum.get_size_dict(code.ixs, tensors)
        optcode = optimize_code(code, size_dict, TreeSA())
        set_logging(true, "data/profile/process_nnz/n=$n.txt")
        value = optcode(tensors...)
    end

    for n in 4:13
        filename = "data/profile/process_nnz/n=$n.txt"
        max_value = find_max_in_file(filename)
        open("data/profile/nnz.txt", "a") do io
            println(io, "n:$n", ", SC_nnz:",log2(max_value))
        end
    end
end


@testset "test num of non-zero elements in the contraction process for sparse tensor network when fixed a config(given by Ip solver) in the encirclement" begin
    for n in [20,22,24]
        L = Int(n/2 - 7)
        t9_lattice_total = generate_TensorNQ_lattice(n)
        pos1, pos0 = fix_one_kind_encirclement(n, L)
        code, dense = generate_masked_3_tensor_network(n,t9_lattice_total,pos1,pos0, Int)
        t1 = @elapsed tensors = map(t -> SparseTensor(t), dense)
        t2 = @elapsed size_dict = OMEinsum.get_size_dict(code.ixs, tensors)
        t3 = @elapsed optcode = optimize_code(code, size_dict, TreeSA())
        cc = OMEinsumContractionOrders.contraction_complexity(optcode, size_dict)
        set_logging(true, "data/profile/process_nnz/sliced_(n-2L)=14_n=$n.txt")
        contraction_stats = @timed value = optcode(tensors...)
        contraction_time_taken = contraction_stats.time  # seconds
        contraction_memory_bytes = contraction_stats.bytes  # bytes
        contraction_memory_mb = contraction_memory_bytes / (1024 * 1024)
        sc = log2(contraction_memory_mb) + 17.0
        
        open("data/profile/sliced_(n-2L)=14.txt", "a") do io
            println(io, "n:$n, L:$L")
            println(io, "Time for sparse map: $(round(t1, digits=6)) seconds")
            println(io, "Time for size_dict: $(round(t2, digits=6)) seconds")
            println(io, "Time for optimize_code: $(round(t3, digits=6)) seconds")
            println(io, "Time for contraction: $(round(contraction_time_taken, digits=6)) seconds")
            println(io, "Memory for contraction: $(round(contraction_memory_mb, digits=6)) MB")
            println(io, "SC: $(round(sc, digits=6))")
            println(io, "contraction complexity: ")
            println(io, cc)
            max_value = find_max_in_file("data/profile/process_nnz/sliced_(n-2L)=14_n=$n.txt")
            println(io, "SC_nnz:",log2(max_value))
            println(io, "--------------------------------")
        end
    end
end