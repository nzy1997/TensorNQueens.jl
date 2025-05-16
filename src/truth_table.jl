using Graphs, OMEinsum, AbstractTrees, TreeWidthSolver
using OMEinsum: getixsv, getiyv, LeafString, uniformsize
using OMEinsum.OMEinsumContractionOrders: IncidenceList, parse_eincode, eo2ct, ContractionTree

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