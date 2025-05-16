using TensorNQueens
using TensorNQueens: generate_tensor, TensorNQ, generate_tensor_network, generate_TensorNQ_lattice,generate_8_tensor_network,generate_3_tensor_network,generate_masked_3_tensor_network, generate_masked_3_bonds, ScNeighborSelector, ScRectangleSelector
using OMEinsum
using OptimalBranching

# Global logging configuration
const EINSUM_LOG_CONFIG = Dict{Symbol, Any}(
    :enabled => false,
    :output => nothing  # can be nothing (for stdout) or a file path
)

"""
    set_logging(enabled::Bool, output_file::Union{String,Nothing}=nothing)

Control whether corresponding logs are written and where the logs are written.
If `output_file` is `nothing`, logs are written to stdout.
"""
function set_logging(enabled::Bool, output_file::Union{String,Nothing}=nothing)
    EINSUM_LOG_CONFIG[:enabled] = enabled
    EINSUM_LOG_CONFIG[:output] = output_file
end

"""
    log(msg::String)

Log a message according to the current logs.
"""
function log(msg::String)
    if EINSUM_LOG_CONFIG[:enabled]
        if EINSUM_LOG_CONFIG[:output] === nothing
            println(msg)
        else
            open(EINSUM_LOG_CONFIG[:output], "a") do io
                println(io, msg)
            end
        end
    end
end


#Each time select a region from pos_vertices=[lattice[i,j].labels[5]], calculate the tbl derived from their mutual constraints and find optimal branching clauses
function position_branching(n::Int, t9_lattice::TensorNQLattice, code::DynamicEinCode, tensors::Vector, sc_target::Int, region_vertices::Vector, solver, scorer)
    #generate the truth table
    region_vertices = [Int(x) for x in region_vertices]
    configs = truth_table(code, tensors, region_vertices)
    if length(configs) == 0
        return [], [], []
    end
    tbl = OptimalBranchingMIS.OptimalBranchingCore.BranchingTable(length(region_vertices),configs)

    #optimal branching
    candidates = OptimalBranchingMIS.OptimalBranchingCore.candidate_clauses(tbl)
    subsets = [OptimalBranchingMIS.OptimalBranchingCore.covered_items(tbl.table, c) for c in candidates]
    num_items = length(tbl.table)
    weights = sc_score_weights(n,t9_lattice,sc_target,region_vertices,candidates,scorer)
    if scorer.delta
        cover = OptimalBranchingMIS.OptimalBranchingCore.weighted_minimum_signed_exact_cover_fixed_point(solver, weights, subsets, num_items, 10.0)
    else
        cover = OptimalBranchingMIS.OptimalBranchingCore.weighted_minimum_signed_exact_cover(solver, weights, subsets, num_items, 10.0)
    end
    
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
#uses ScRectangleSelector, which selects regions based on the sc reduction
#branching untill the sc is less than sc_target
function tensor_branching(n::Int, t9_lattice::TensorNQLattice, pos1::Vector, pos0::Vector, coefficient::Float64, sc_target::Int, region_selector::ScRectangleSelector, IP_solver, scorer, depth::Int)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    cc = contraction_complexity(optcode, uniformsize(optcode, 2))
    if cc.sc <= sc_target
        counting_branch = optcode(tensors...)[]
        # counting_branch = 1   # When we only want to know the number of slices and tc, we don't actually perform tn contraction
        log("\nend:")
        log("\tdepth: $depth, coefficient: $coefficient")
        log("\tconfig: pos1: $(unique(pos1)), pos0: $(unique(pos0))")
        log("\tsc: $(cc.sc)")
        log("\tcounting_branch: $counting_branch")
        return [cc], [coefficient*counting_branch]
    end
    ccs = []
    countings = []
    region_vertices = branching_region(n,t9_lattice,region_selector,code,optcode)
    branches, branch_coefficients, branch_weights = position_branching(n, t9_lattice, code, tensors, sc_target, region_vertices, IP_solver, scorer)
    log("\ndepth: $depth, coefficient: $coefficient")
    log("branches: $branches")
    log("sc before branch: $(cc.sc)")
    for branch_id in 1:length(branches)
        branch_coefficient = branch_coefficients[branch_id] * coefficient
        branch_pos1 = vcat(pos1, branches[branch_id][1])
        branch_pos0 = vcat(pos0, branches[branch_id][2])
        ccs_branch, countings_branch = tensor_branching(n, t9_lattice, branch_pos1, branch_pos0, branch_coefficient, sc_target, region_selector, IP_solver, scorer, depth+1)
        append!(ccs, ccs_branch)
        append!(countings, countings_branch)
    end
    if depth == 1
        log("\nterminate: ")
        log("ccs: ")
        for cc in ccs
            log("$cc")
            log("_______________________________________")
        end
        log("countings: $countings")
        log("total tc: $(log2(sum([2^(cc.tc) for cc in ccs])))")
        log("countings sum: $(sum(countings))")
        log("slice number: $(length(countings))")
    end
    return ccs, countings
end


#uses ScRectangleD4Selector, which selects regions based on the reduction of the number of rank>4 tensors
#branching untill the distance is less than some given constant
function tensor_branching(n::Int, t9_lattice::TensorNQLattice, pos1::Vector, pos0::Vector, coefficient::Float64, sc_target::Int, region_selector::ScRectangleD4Selector, IP_solver, scorer, depth::Int)
    pos_vertices = vec([t9_lattice.lattice[i,j].labels[5] for i in 1:n, j in 1:n])
    distance = distance_square_lattice(n, t9_lattice, pos_vertices, pos1, pos0)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    cc = contraction_complexity(optcode, uniformsize(optcode, 2))
    if cc.sc <= sc_target
        counting_branch = optcode(tensors...)[]
        # counting_branch = 1
        log("\nend:")
        log("\tdepth: $depth, coefficient: $coefficient")
        log("\tconfig: pos1: $(unique(pos1)), pos0: $(unique(pos0))")
        log("\tcounting_branch: $counting_branch")
        log("\tdistance: $distance")
        log("cc: \n$cc")
        return [distance], [cc.sc], [coefficient*counting_branch]
    end

    distances = []
    countings = []
    scs = []
    region_vertices = branching_region(n,t9_lattice,region_selector,pos1,pos0,code)
    branches, branch_coefficients, branch_weights = position_branching(n, t9_lattice, code, tensors, sc_target, region_vertices, IP_solver, scorer)
    log("\ndepth: $depth, coefficient: $coefficient")
    log("branches: $branches")
    log("distance before branch: $distance")
    log("sc before branch: $(cc.sc)")
    for branch_id in 1:length(branches)
        branch_coefficient = branch_coefficients[branch_id] * coefficient
        branch_pos1 = unique(vcat(pos1, branches[branch_id][1]))
        branch_pos0 = unique(vcat(pos0, branches[branch_id][2]))
        branches_branch, scs_branch, countings_branch = tensor_branching(n, t9_lattice, branch_pos1, branch_pos0, branch_coefficient, sc_target, region_selector, IP_solver, scorer, depth+1)
        append!(distances, branches_branch)
        append!(countings, countings_branch)
        append!(scs, scs_branch)
    end
    if depth == 1
        log("\nterminate: ")
        log("distances: $distances")
        log("countings: $countings")
        log("scs: $scs")
        log("countings sum: $(sum(countings))")
        log("slice number: $(length(countings))")
    end
    return distances, scs, countings
end


#In every branching step, choose the row with the least survival length l, and do the naive l-branches branching
function naive_tensor_branching(n::Int, t9_lattice::TensorNQLattice, pos1::Vector, pos0::Vector, sc_target::Int, depth::Int)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    cc = contraction_complexity(optcode, uniformsize(optcode, 2))
    if cc.sc <= sc_target
        counting_branch = optcode(tensors...)[]
        # counting_branch = 1
        log("\nend:")
        log("\tdepth: $depth")
        log("\tsc: $(cc.sc)")
        log("\tcounting_branch: $counting_branch")
        return [cc], [counting_branch]
    end
    ccs = []
    countings = []
    bonds = vcat(getixsv(code)...)
    min_survival_length = Inf
    chosen_i = 0
    chosen_j = 0
    for i in 1:n
        survival_length = count(j -> t9_lattice.lattice[i,j].labels[5] in bonds, 1:n)
        if survival_length < min_survival_length && survival_length > 0
            min_survival_length = survival_length
            chosen_i = i
        end
    end
    for j in 1:n
        survival_length = count(i -> t9_lattice.lattice[i,j].labels[5] in bonds, 1:n)
        if survival_length < min_survival_length && survival_length > 0
            min_survival_length = survival_length
            chosen_j = j
            chosen_i = 0
        end
    end
    branches = []
    if chosen_i != 0
        js = [j for j in 1:n if t9_lattice.lattice[chosen_i,j].labels[5] in bonds]
        for j in js
            push!(branches, [[(chosen_i,j)], [(chosen_i,jj) for jj in js if jj != j]])
        end
    else
        is = [i for i in 1:n if t9_lattice.lattice[i,chosen_j].labels[5] in bonds]
        for i in is
            push!(branches, [[(i,chosen_j)], [(ii,chosen_j) for ii in is if ii != i]])
        end
    end
    log("\ndepth: $depth")
    log("branches: $branches")
    log("chosen_i: $chosen_i, chosen_j: $chosen_j")
    log("sc before branch: $(cc.sc)")
    for branch_id in 1:length(branches)
        branch_pos1 = vcat(pos1, branches[branch_id][1])
        branch_pos0 = vcat(pos0, branches[branch_id][2])
        ccs_branch, countings_branch = naive_tensor_branching(n, t9_lattice, branch_pos1, branch_pos0, sc_target, depth+1)
        append!(ccs, ccs_branch)
        append!(countings, countings_branch)
    end
    if depth == 1
        log("\nterminate: ")
        log("ccs: ")
        for cc in ccs
            log("$cc")
            log("_______________________________________")
        end
        log("countings: $countings")
        log("total tc: $(log2(sum([2^(cc.tc) for cc in ccs])))")
        log("countings sum: $(sum(countings))")
        log("slice number: $(length(countings))")
    end
    return ccs, countings
end


#In every branching step, choose a random survival position, and do the naive 2-branches(including i, excluding i) branching
function random_naive_tensor_branching(n::Int, t9_lattice::TensorNQLattice, pos1::Vector, pos0::Vector, sc_target::Int, depth::Int)
    code, tensors = generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,Int)
    optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
    cc = contraction_complexity(optcode, uniformsize(optcode, 2))
    if cc.sc <= sc_target
        counting_branch = optcode(tensors...)[]
        # counting_branch = 1
        log("\nend:")
        log("\tdepth: $depth")
        log("\tsc: $(cc.sc)")
        log("\tcounting_branch: $counting_branch")
        return [cc], [counting_branch]
    end
    ccs = []
    countings = []
    bonds = vcat(getixsv(code)...)
    survival_poses = [(i,j) for i in 1:n for j in 1:n if t9_lattice.lattice[i,j].labels[5] in bonds]
    chosen_pos = rand(survival_poses)
    branches = [[[chosen_pos],[]], [[],[chosen_pos]]]
    log("\ndepth: $depth")
    log("branches: $branches")
    log("chosen_pos: $chosen_pos")
    log("sc before branch: $(cc.sc)")
    for branch_id in 1:length(branches)
        branch_pos1 = vcat(pos1, branches[branch_id][1])
        branch_pos0 = vcat(pos0, branches[branch_id][2])
        ccs_branch, countings_branch = random_naive_tensor_branching(n, t9_lattice, branch_pos1, branch_pos0, sc_target, depth+1)
        append!(ccs, ccs_branch)
        append!(countings, countings_branch)
    end
    if depth == 1
        log("\nterminate: ")
        log("ccs: ")
        for cc in ccs
            log("$cc")
            log("_______________________________________")
        end
        log("countings: $countings")
        log("total tc: $(log2(sum([2^(cc.tc) for cc in ccs])))")
        log("countings sum: $(sum(countings))")
        log("slice number: $(length(countings))")
    end
    return ccs, countings
end