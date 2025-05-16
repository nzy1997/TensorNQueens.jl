using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, generate_masked_3_bonds, ScNeighborSelector, ScRectangleSelector
using OptimalBranching
using Graphs, OMEinsum, AbstractTrees, TreeWidthSolver
using OMEinsum: getixsv, getiyv, LeafString, uniformsize
using OMEinsum.OMEinsumContractionOrders: IncidenceList, parse_eincode, eo2ct, ContractionTree

struct ExactScScorer
    sc_target::Int  #target sc
end

struct TruncateBondScorer
    bond_limit::Int   #maximum rank of the tensors
    weighted::Bool  #whether to weight the score where the contribution = 2^(bond_rank - bond_limit)
    delta::Bool  #whether to use delta 
end


function _log2_einsize(eincode::ET, size_dict::Dict{LT, Int}) where {ET, LT}
    return foldl((x, y) -> x + log2(size_dict[y]), eincode.iy, init = 0.0)
end


function list_subtree(code::CT, size_dict::Dict{LT, Int}, threshold::T) where {CT, LT, T}
    subtrees = Vector{CT}()
    for subtree in PostOrderDFS(code)
        (subtree isa LeafString) && continue
        if _log2_einsize(subtree.eins, size_dict) ≥ threshold
            push!(subtrees, subtree)
        end
    end
    return subtrees
end

#Calculate the weight of each candidate, which is the truncated_sc of the new TN after fixing the corresponding variable values of the candidate
function sc_score_weight(n::Int, t9_lattice::TensorNQLattice, sc_target::Int, region_vertices::Vector{Int}, candidate, scorer::ExactScScorer)
    pos0, pos1 = generate_pos_vec_local(t9_lattice, region_vertices, candidate.mask, candidate.val)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    large_tensors = list_subtree(optcode.eins, uniformsize(optcode, 2), 2)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    #sc of all tree bags bigger than sc_target
    score = sum(2.0^(max(0, length(lt_iy) - sc_target)) - 1.0 for lt_iy in large_tensors_iys) + 1
    return score
end


function sc_score_weights(n::Int, t9_lattice::TensorNQLattice, sc_target::Int,  region_vertices::Vector{Int}, candidates::Vector, scorer::ExactScScorer)
    weights = [sc_score_weight(n, t9_lattice, sc_target, region_vertices, c, scorer) for c in candidates]
    return weights
end


#score based on the distance to the square lattice tn (with all tensor ranks equal to 4)
function sc_score_weights(n::Int, t9_lattice::TensorNQLattice, sc_target::Int, region_vertices::Vector{Int}, candidates::Vector, scorer::TruncateBondScorer)
    pos_vertices =vec([t9_lattice.lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    weights = []
    bonds_original = generate_masked_3_bonds(n, t9_lattice, [], [], Int)
    rank_original = sum([length(unique(Iterators.flatten(bond for bond in bonds_original if pos_vertex in bond))) - scorer.bond_limit for pos_vertex in pos_vertices])
            
    for candidate in candidates
        pos0, pos1 = generate_pos_vec_local(t9_lattice, region_vertices, candidate.mask, candidate.val)
        bonds_candidate = generate_masked_3_bonds(n, t9_lattice, pos1, pos0, Int)
        weight = 0.0
        for pos_vertex in pos_vertices
            rank = length(unique(Iterators.flatten(bond for bond in bonds_candidate if pos_vertex in bond)))
            if scorer.weighted
                weight += max(1,2.0^(rank - scorer.bond_limit)) - 1
            else
                weight += max(0,rank - scorer.bond_limit)
            end
        end
        if scorer.delta
            push!(weights, rank_original-weight)
        else
            push!(weights, weight)
        end
    end
    return weights
end