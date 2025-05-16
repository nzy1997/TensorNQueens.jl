using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, generate_masked_3_bonds, ScNeighborSelector, ScRectangleSelector
using OMEinsum
using OptimalBranching

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

struct ScRectangleD4Selector 
    k_ud::Int   #k-layers of neighbors in the up-down direction
    k_lr::Int  #k-layers of neighbors in the left-right direction
    n_max::Int  #maximum number of vertices in the region
    sc_target::Int 
end

#When region_selector == ScNeighborSelector:
#Generate a graph for pos_vertices based on their connections in the code, and generate k-neighbors for each vertex
function generate_neighbors(n::Int, t9_lattice::TensorNQLattice, pos_vertices::Vector{Int}, code::DynamicEinCode, region_selector::ScNeighborSelector)
    lattice = t9_lattice.lattice
    bonds = getixsv(code)
    pos_vertices_neighbors = [
        unique(Iterators.flatten(bond for bond in bonds if pos_vertex in bond))
        for pos_vertex in pos_vertices
    ]
    graph = SimpleGraph(length(pos_vertices))

    for i in 1:length(pos_vertices)
        for j in i+1:length(pos_vertices)
            if !isempty(intersect(pos_vertices_neighbors[i], pos_vertices_neighbors[j]))
                add_edge!(graph,i,j)
            end
        end
    end
    neighbors = []
    for v in 1:length(pos_vertices)
        neighbor = OptimalBranchingMIS.neighbor_cover(graph,v,region_selector.k)[1]
        neighbor = [pos_vertices[i] for i in neighbor]
        push!(neighbors,neighbor)
    end
    return neighbors
end

#When region_selector == ScRectangleSelector:
#Generate the rectangle region of size [k_ud,k_lr] centered at each pos_vertex
function generate_neighbors(n::Int, t9_lattice::TensorNQLattice, pos_vertices::Vector{Int}, code::DynamicEinCode, region_selector::Union{ScRectangleSelector, ScRectangleD4Selector})
    neighbors = []
    for pos_vertex in pos_vertices
        iv,jv = Tuple(findfirst(x -> x.labels[5] == pos_vertex, t9_lattice.lattice))
        neighbor = []
        for i in max(1,iv-region_selector.k_ud):min(n,iv+region_selector.k_ud)
            for j in max(1,jv-region_selector.k_lr):min(n,jv+region_selector.k_lr)
                if t9_lattice.lattice[i,j].labels[5] in pos_vertices
                    push!(neighbor,t9_lattice.lattice[i,j].labels[5])
                end
            end
        end
        push!(neighbors,neighbor)
    end
    return neighbors
end

#When region_selector == ScNeighborSelector:
#(1) Generate a graph for pos_vertices based on their connections in the code
#(2) Find k-th order neighbors for each pos_vertex on this graph
#(3) On the contraction tree corresponding to optcode, remove vertices from each neighbor and see which gives the smallest truncated_sc
function branching_region(n::Int, t9_lattice::TensorNQLattice, region_selector::ScNeighborSelector, code, optcode)
    pos_vertices = vec([t9_lattice.lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    bonds = getixsv(code)
    vertices_in_bonds = unique(Iterators.flatten(bonds))
    pos_vertices = filter(v -> v in vertices_in_bonds, pos_vertices)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)
    large_tensors = list_subtree(optcode.eins, uniformsize(optcode, 2), 2)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    min_score = Inf
    min_neighbor = [] 
    for neighbor in neighbors
        filtered_tensors_iys = [setdiff(lt_iy, neighbor) for lt_iy in large_tensors_iys]
        score = sum(2.0^(max(0, length(lt_iy) - region_selector.sc_target)) - 1.0 for lt_iy in filtered_tensors_iys) + 1
        if score < min_score
            min_score = score
            min_neighbor = neighbor
        end
    end
    return min_neighbor
end
    

#When region_selector == ScRectangleSelector:
#Generate the rectangle region of size [k_ud,k_lr] centered at each pos_vertex
function branching_region(n::Int, t9_lattice::TensorNQLattice, region_selector::ScRectangleSelector, code, optcode)
    pos_vertices = vec([t9_lattice.lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    bonds = getixsv(code)
    vertices_in_bonds = unique(Iterators.flatten(bonds))
    pos_vertices = filter(v -> v in vertices_in_bonds, pos_vertices)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)
    large_tensors = list_subtree(optcode.eins, uniformsize(optcode, 2), 2)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    min_score = Inf
    min_neighbor = [] 
    for neighbor in neighbors
        filtered_tensors_iys = [setdiff(lt_iy, neighbor) for lt_iy in large_tensors_iys]
        score = sum(2.0^(max(0, length(lt_iy) - region_selector.sc_target)) - 1.0 for lt_iy in filtered_tensors_iys) + 1
        if score < min_score
            min_score = score
            min_neighbor = neighbor
        end
    end
    return min_neighbor
end


function distance_square_lattice(n::Int, t9_lattice::TensorNQLattice, pos_vertices::Vector, pos1::Vector, pos0::Vector)
    bonds = generate_masked_3_bonds(n, t9_lattice, pos1, pos0, Int)
    L1_distance = sum([max(0,length(unique(Iterators.flatten(bond for bond in bonds if pos_vertex in bond))) - 4) for pos_vertex in pos_vertices])
    return L1_distance
end


#When region_selector == ScRectangleD4Selector:
#Generate the rectangle region of size [k_ud,k_lr] centered at each pos_vertex
#The score of each rectangle is based on the distance of the current tn to square lattice tn
function branching_region(n::Int, t9_lattice::TensorNQLattice, region_selector::ScRectangleD4Selector, pos1::Vector, pos0::Vector, code)
    pos_vertices = vec([t9_lattice.lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    bonds = getixsv(code)
    vertices_in_bonds = unique(Iterators.flatten(bonds))
    pos_vertices = filter(v -> v in vertices_in_bonds, pos_vertices)
    neighbors = generate_neighbors(n,t9_lattice,pos_vertices,code,region_selector)

    min_score = Inf
    min_neighbor = [] 
    pos1_original = copy(pos1)
    pos0_original = copy(pos0)
    for i in eachindex(neighbors)
        neighbor = neighbors[i]
        pos1 = copy(pos1_original)
        pos0 = copy(pos0_original)
        append!(pos1, [Tuple(findfirst(x -> x.labels[5] == pos_vertices[i], t9_lattice.lattice))])
        append!(pos0, [Tuple(findfirst(x -> x.labels[5] == j, t9_lattice.lattice)) for j in neighbor if j != pos_vertices[i]])
        pos1 = unique(pos1)
        pos0 = unique(pos0)
        score = distance_square_lattice(n, t9_lattice, pos_vertices, pos1, pos0)
        if score < min_score
            min_score = score
            min_neighbor = neighbor
        end
    end
    return min_neighbor
end
