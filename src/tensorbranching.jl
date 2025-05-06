using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, ScNeighborSelector, ScRectangleSelector
using OMEinsum
using OptimalBranching


function truth_table(code::DynamicEinCode, tensors::Vector, region_vertices::Vector{Int})
    bonds = getixsv(code)
    involved_bonds = Int[]   #bonds involved in the local tn contraction
    local_bonds = Vector{Int}[]  #final bonds in the local tn contraction
    local_tensors = []  #final tensors in the local tn contraction

    #For 3-order tensors, only those containing region_vertices will be included in local tn contraction
    for tid in eachindex(bonds)
        bond = bonds[tid]
        tensor = tensors[tid]
        if length(bond) == 3 && any(v -> v in bond, region_vertices)
            append!(involved_bonds, bond)
            push!(local_bonds, bond)
            push!(local_tensors, tensor)
        end
    end
    involved_bonds = unique(involved_bonds) 

    #For 1-order tensors, only those containing bonds involved in the 3-order tensors above will be included in local tn contraction
    for tid in eachindex(bonds)
        bond = bonds[tid]
        tensor = tensors[tid]
        if length(bond) == 1 && bond[1] in involved_bonds && !(bond[1] in region_vertices)
            push!(local_bonds, bond)
            push!(local_tensors, tensor)
            deleteat!(involved_bonds, findfirst(==(bond[1]), involved_bonds))
        end
    end

    # Only constraints between region_vertices take effect here. For constraints from outside, except for those on the boundary of the whole model,
    # all directions from outside the region are free, so we only need to connect a (1,1)
    for bond in involved_bonds
        push!(local_bonds, [bond])
        push!(local_tensors, [one(Int),one(Int)])
    end

    # Add an open leg on all region_vertices
    code = DynamicEinCode(local_bonds, region_vertices)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    alpha_tensor = optcode(local_tensors...)
    configs_cartesian = findall(x->x!=0, alpha_tensor)
    configs = [[collect(Tuple(c)) .- 1] for c in configs_cartesian]
    return configs
end

#Calculate the weight of each candidate, which is the truncated_sc of the new TN after fixing the corresponding variable values of the candidate
function sc_score_weight(n::Int, t9_lattice::TensorNQLattice, sc_target::Int, region_vertices::Vector{Int}, candidate)
    pos0, pos1 = generate_pos_vec_local(t9_lattice, region_vertices, candidate.mask, candidate.val)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    large_tensors = list_subtree(optcode.eins, uniformsize(optcode, 2), 2)
    large_tensors_iys = [Set(t.eins.iy) for t in large_tensors]

    #sc of all tree bags bigger than sc_target
    score = sum(2.0^(max(0, length(lt_iy) - sc_target)) - 1.0 for lt_iy in large_tensors_iys) + 1
    return score
end


function sc_score_weights(n::Int, t9_lattice::TensorNQLattice, sc_target::Int,  region_vertices::Vector{Int}, candidates)
    weights = [sc_score_weight(n, t9_lattice, sc_target, region_vertices, c) for c in candidates]
    return weights
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
function generate_neighbors(n::Int, t9_lattice::TensorNQLattice, pos_vertices::Vector{Int}, code::DynamicEinCode, region_selector::ScRectangleSelector)
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
    pos_vertices =vec([t9_lattice.lattice[i,j].labels[5] for i in 1:n, j in 1:n])
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


#Each time select a region from pos_vertices=[lattice[i,j].labels[5]], calculate the tbl derived from their mutual constraints and find optimal branching clauses
function position_branching(n::Int, t9_lattice::TensorNQLattice, code::DynamicEinCode, tensors::Vector, sc_target::Int, region_vertices::Vector, solver)
    #generate the truth table
    region_vertices = [Int(x) for x in region_vertices]
    configs = truth_table(code, tensors, region_vertices)
    tbl = OptimalBranchingMIS.OptimalBranchingCore.BranchingTable(length(region_vertices),configs)

    #optimal branching
    candidates = OptimalBranchingMIS.OptimalBranchingCore.candidate_clauses(tbl)
    subsets = [OptimalBranchingMIS.OptimalBranchingCore.covered_items(tbl.table, c) for c in candidates]
    num_items = length(tbl.table)
    weights = sc_score_weights(n,t9_lattice,sc_target,region_vertices,candidates)
    cover = OptimalBranchingMIS.OptimalBranchingCore.weighted_minimum_signed_exact_cover(solver, weights, subsets, num_items, 10.0)
    
    #pick the valid branches
    picked_scs = findall(x -> (x > 1e-6) || (x < -1e-6), cover)
    branch_coefficients = cover[picked_scs] #coefficients of the branches contributed to the whole counting
    clauses = candidates[picked_scs]
    branch_weights = weights[picked_scs]
    branches = []
    for i in 1:length(clauses)
        pos0_branch, pos1_branch = generate_pos_vec_local(t9_lattice, region_vertices, clauses[i].mask, clauses[i].val)
        push!(branches, [pos1_branch,pos0_branch])
    end
    return branches, branch_coefficients, branch_weights
end


#Depth-first search, where branch_coefficients are the coefficients to multiply for each branch
function tensor_branching(n::Int, t9_lattice::TensorNQLattice, pos1::Vector, pos0::Vector, coefficient::Float64, sc_target::Int, region_selector::ScRectangleSelector, IP_solver)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    cc = contraction_complexity(optcode, uniformsize(optcode, 2))
    if cc.sc <= sc_target
        counting_branch = optcode(tensors...)[]
        return [cc], [coefficient*counting_branch]
    end
    ccs = []
    countings = []
    region_vertices = branching_region(n,t9_lattice,region_selector,code,optcode)
    branches, branch_coefficients, branch_weights = position_branching(n, t9_lattice, code, tensors, sc_target, region_vertices, IP_solver)
    println("branches:", branches)
    println("cc before branch:", cc)
    for branch_id in 1:length(branches)
        branch_coefficient = branch_coefficients[branch_id] * coefficient
        branch_pos1 = vcat(pos1, branches[branch_id][1])
        branch_pos0 = vcat(pos0, branches[branch_id][2])
        ccs_branch, countings_branch = tensor_branching(n, t9_lattice, branch_pos1, branch_pos0, branch_coefficient, sc_target, region_selector, IP_solver)
        append!(ccs, ccs_branch)
        append!(countings, countings_branch)
    end
    return ccs, countings
end