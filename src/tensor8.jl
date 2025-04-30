# 3 4 7
# 1   8
# 2 5 6
function generate_8_tensor(T)
    t = generate_tensor(T)
    t10,t01,t11 = generate01tensors(T)
    return ein"abcdefghi,e -> abcdfghi"(t,t11)
end

function get_pos8()
    return [(0,-1), (1,-1), (-1,-1), (-1,0), (1,0), (1,1),(-1,1),(0,1)]
end

function generate_8_TensorNQ_lattice(n::Int)
    # res = Matrix{TensorNQ}(undef,n,n)
    res = fill(TensorNQ(fill(-1,9)),n,n)
    label_count = 0

    pos_vec = get_pos8()

    pos11 = Int[]
    pos01 = Int[]
    pos10 = Int[]
    for i in 1:n
        for j in 1:n
            index_vec = fill(-1,8)
            for (ind,pos) in enumerate(pos_vec)
                if checkbounds(Bool,res,i+pos[1],j+pos[2])
                    if res[i+pos[1],j+pos[2]].labels[9-ind] != -1
                        index_vec[ind] = res[i+pos[1],j+pos[2]].labels[9-ind]
                    else
                        label_count += 1
                        index_vec[ind] = label_count
                    end
                else
                    label_count += 1
                    index_vec[ind] = label_count
                    if ind < 5
                        push!(pos10,label_count)
                    elseif ind == 8 || ind == 5
                        push!(pos01,label_count)
                    else
                        push!(pos11,label_count)
                    end
                end
            end
            res[i,j] = TensorNQ(index_vec)
        end
    end
    return res, pos10,pos01,pos11,label_count
end

function generate_8_tensor_network(n::Int,T)
    lattice,pos10,pos01,pos11,_ = generate_8_TensorNQ_lattice(n)
    t8_ixs = getfield.(vec(lattice),:labels)
    t = generate_8_tensor(T)
    t10,t01,t11 = generate01tensors(T)
    return DynamicEinCode(t8_ixs ∪ [[p] for p in pos10] ∪ [[p] for p in pos01] ∪  [[p] for p in pos11],Int[]),[fill(t,length(t8_ixs))...,fill(t10,length(pos10))...,fill(t01,length(pos01))...,fill(t11,length(pos11))...]
end