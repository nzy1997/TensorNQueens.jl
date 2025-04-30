function generate_3_tensor(T)
    return T[[1 0;0 1];;;[0 1;0 0]]
end


function generate_3_tensor_network(n::Int,T)
    t9_lattice = generate_TensorNQ_lattice(n)
    return generate_3_tensor_network(t9_lattice,T)
end
function generate_3_tensor_network(t9_lattice::TensorNQ_lattice,T)
    lattice,pos10,pos01,pos11 = t9_lattice.lattice,t9_lattice.pos10,t9_lattice.pos01,t9_lattice.pos11
    t3_ixs = Vector{Vector{Int}}()
    t9_ixs = getfield.(vec(lattice),:labels)

    for vec in t9_ixs
        for i in 1:4
            push!(t3_ixs,[vec[i],vec[10-i],vec[5]])
        end
    end
    t10,t01,t11 = generate01tensors(T)
    t3 = generate_3_tensor(T)
    return DynamicEinCode(t3_ixs ∪ [[p] for p in pos10] ∪ [[p] for p in pos01] ∪  [[p] for p in pos11],Int[]),[fill(t3,length(t3_ixs))...,fill(t10,length(pos10))...,fill(t01,length(pos01))...,fill(t11,length(pos11))...]
end


function generate_masked_3_tensor_network(t9_lattice::TensorNQ_lattice,pos1::Vector{Tuple{Int,Int}},pos0::Vector{Tuple{Int,Int}},T)
    lattice,pos10,pos01,pos11 = t9_lattice.lattice,t9_lattice.pos10,t9_lattice.pos01,t9_lattice.pos11
    pos10 = copy(pos10)
    pos01 = copy(pos01)
    pos11 = copy(pos11)
    lattice_copy = deepcopy(lattice)
    for (i,j) in pos1
        for (index,(id,jd)) in enumerate(get_pos8())
            i_new,j_new = i,j
            while checkbounds(Bool,lattice,i_new+id,j_new+jd)
                i_new += id
                j_new += jd
                push!(pos0,(i_new,j_new))
            end
            rm_t = lattice[i_new,j_new].labels[index]
            if index < 5
                setdiff!(pos10,rm_t)
            elseif index == 8 || index == 5
                setdiff!(pos01,rm_t)
            else
                setdiff!(pos11,rm_t)
            end
        end
    end

    for (i,j) in pos0
        for (index,(id,jd)) in enumerate(get_pos8())
            i_new,j_new = i,j
            while checkbounds(Bool,lattice,i_new+id,j_new+jd)
                i_new += id
                j_new += jd
            end
        end
    end
    t3_ixs = Vector{Vector{Int}}()
    t9_ixs = getfield.(vec(lattice),:labels)

    for vec in t9_ixs
        for i in 1:4
            push!(t3_ixs,[vec[i],vec[10-i],vec[5]])
        end
    end
    t10,t01,t11 = generate01tensors(T)
    t3 = generate_3_tensor(T)
    return DynamicEinCode(t3_ixs ∪ [[p] for p in pos10] ∪ [[p] for p in pos01] ∪  [[p] for p in pos11],Int[]),[fill(t3,length(t3_ixs))...,fill(t10,length(pos10))...,fill(t01,length(pos01))...,fill(t11,length(pos11))...]
end