function generate_3_tensor(T)
    return T[[1 0;0 1];;;[0 1;0 0]]
end

function generate_3_tensor_network(n::Int,T)
    t9_lattice = generate_TensorNQ_lattice(n)
    return generate_3_tensor_network(t9_lattice,T)
end
function generate_3_tensor_network(t9_lattice::TensorNQLattice,T)
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

function generate_masked_3_tensor_network(n,t9_lattice::TensorNQLattice,pos1::Vector,pos0::Vector,T)
    masked_lattice = generate_MaskedTensorNQLattice(n,t9_lattice,pos1,pos0,T)
    t3_ixs = Vector{Vector{Int}}()

    for i in 1:n
        for j in 1:n
            if (i,j) ∉ masked_lattice.pos0 && (i,j) ∉ masked_lattice.pos1
                for index in 1:4
                    push!(t3_ixs,[masked_lattice.lattice[i,j].labels[index],masked_lattice.lattice[i,j].labels[10-index],masked_lattice.lattice[i,j].labels[5]])
                end
            end
        end
    end
    t10,t01,t11 = generate01tensors(T)
    t3 = generate_3_tensor(T)
    # show_lattice(lattice)
    # show_lattice(lattice_copy)
    # @show pos10
    # @show pos01
    # @show pos11
    return DynamicEinCode(vcat(t3_ixs , [[p] for p in masked_lattice.pos10] , [[p] for p in masked_lattice.pos01] ,  [[p] for p in masked_lattice.pos11]),Int[]),[fill(t3,length(t3_ixs))...,fill(t10,length(masked_lattice.pos10))...,fill(t01,length(masked_lattice.pos01))...,fill(t11,length(masked_lattice.pos11))...]
end

function show_lattice(lattice::Matrix{TensorNQ})
    n = size(lattice,1)
    mat = zeros(Int,3*n,3*n)
    for i in 1:n
        for j in 1:n
            mat[3*i-2:3*i,3*j-2:3*j] = [lattice[i,j].labels[3] lattice[i,j].labels[4] lattice[i,j].labels[8]; lattice[i,j].labels[1] lattice[i,j].labels[5] lattice[i,j].labels[9]; lattice[i,j].labels[2] lattice[i,j].labels[6] lattice[i,j].labels[7]]
        end
    end
    display(mat)
end

function remove1!(pos0,pos10,pos01,pos11,lattice_copy,i,j)
    for (index,(id,jd)) in enumerate(get_pos8())
        i_new,j_new = i,j
        while checkbounds(Bool,lattice_copy,i_new+id,j_new+jd)
            i_new += id
            j_new += jd
            push!(pos0,(i_new,j_new))
        end
        index2 = index > 4 ? index + 1 : index
        rm_t = lattice_copy[i_new,j_new].labels[index2]
        if index2 < 5
            setdiff!(pos10,rm_t)
        elseif index2 == 9 || index2 == 6
            setdiff!(pos01,rm_t)
        else
            setdiff!(pos11,rm_t)
        end
    end
    setdiff!(pos11,lattice_copy[i,j].labels[5])
end

function remove0!(pos01,pos11,lattice_copy,i,j)
    for (index,(id,jd)) in enumerate(get_pos8()[1:4])
        i_new,j_new = i+id,j+jd
        if checkbounds(Bool,lattice_copy,i_new,j_new)
            edge_ind = lattice_copy[i_new,j_new].labels[10 - index]
            lattice_copy[i,j].labels[index] = edge_ind
            
        else
            edge_ind = lattice_copy[i,j].labels[index]

        end
        i_new,j_new = i - id,j - jd
        if checkbounds(Bool,lattice_copy,i_new,j_new)
            lattice_copy[i_new,j_new].labels[index] = edge_ind
        else
            if index == 1 || index == 4
                replace!(pos01,lattice_copy[i,j].labels[10-index] => edge_ind)
            else
                replace!(pos11,lattice_copy[i,j].labels[10-index] => edge_ind)
            end
        end
        lattice_copy[i,j].labels[10 - index] = edge_ind
    end
    setdiff!(pos11,lattice_copy[i,j].labels[5])
end

# label the chess board 
# 1 n+1 ... n^2-n +1
# |  |      |    
# 2 n+2 ... n^2-n+2
# |  |      |    
# ... ... ...
# n 2n  ... n^2

function generate_pos_vec(n,mask,val)
    pos0 = Vector{Tuple{Int,Int}}()
    pos1 = Vector{Tuple{Int,Int}}()
    for i in 1:n
        for j in 1:n
            bit_pos = (i-1)*n + j
            if readbit(mask,bit_pos) == 1
                if readbit(val,bit_pos) == 1
                    push!(pos1,(j,i))
                else
                    push!(pos0,(j,i))
                end
            end
        end
    end
    return pos0,pos1
end

function generate_masked_3_tensor_network(n,t9_lattice::TensorNQLattice,mask,val,T)
    pos0,pos1 = generate_pos_vec(n,mask,val)
    return generate_masked_3_tensor_network(n,t9_lattice,pos1,pos0,T)
end