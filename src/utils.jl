using Graphs, OMEinsum, AbstractTrees, TreeWidthSolver
using OMEinsum: getixsv, getiyv, LeafString, uniformsize
using OMEinsum.OMEinsumContractionOrders: IncidenceList, parse_eincode, eo2ct, ContractionTree

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


struct ScNeighborSelector 
    k::Int  #k-layers of neighbors
    n_max::Int  #maximum number of vertices in the region
    sc_target::Int 
end


struct ScRectangleSelector 
    k_ud::Int   #k-layers of neighbors in the up-down direction
    k_lr::Int  #k-layers of neighbors in the left-right direction
    n_max::Int  #maximum number of vertices in the region
    sc_target::Int 
end