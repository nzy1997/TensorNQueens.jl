function generate_planar_lattice(n::Int)
    res = fill(TensorNQ(fill(-1,8)),n,n)
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
                    if res[i+pos[1],j+pos[2]].labels[9-ind] != -1 && (ind == 1 || ind == 4 || ind == 5 || ind == 8)
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
    return TensorNQLattice(res, pos10,pos01,pos11),label_count
end

# 1 3
#  X
# 2 4
function swap_tensor(T)
    res = zeros(T, fill(2,4)...)
    res[1,1,1,1] = one(T)
    res[2,1,1,2] = one(T)
    res[1,2,2,1] = one(T)
    res[2,2,2,2] = one(T)
    return res
end

function generate_planar_tensor_network(n::Int,T)
    tl,max_label = generate_planar_lattice(n)
    t = swap_tensor(T)
    t10,t01,t11 = generate01tensors(T)
    t8 = generate_8_tensor(T)
    tensors = PlanarTensor{Int,T}[]
    for j in 1:n
        push!(tensors,PlanarTensor(t10,[tl.lattice[1,j].labels[4]], Float64(j), 0.0))
        push!(tensors,PlanarTensor(t10,[tl.lattice[1,j].labels[3]], Float64(j)-0.3, 0.0))
        push!(tensors,PlanarTensor(t11,[tl.lattice[1,j].labels[7]], Float64(j)+ 0.3, 0.0))
    end

    for i in 1:n
        push!(tensors,PlanarTensor(t10,[tl.lattice[i,1].labels[1]], 0.0, Float64(i)))
        push!(tensors,PlanarTensor(t10,[tl.lattice[i,1].labels[2]], 0.0, Float64(i) + 0.3))
        if i > 1
            push!(tensors,PlanarTensor(t10,[tl.lattice[i,1].labels[3]], 0.0, Float64(i) - 0.3))
        end
    end

    for i in 1:n
        for j in 1:n
            push!(tensors,PlanarTensor(t8,tl.lattice[i,j].labels, Float64(j), Float64(i)))
            if i < n && j < n
                push!(tensors,PlanarTensor(t,[tl.lattice[i,j].labels[6],tl.lattice[i+1,j].labels[7],tl.lattice[i,j+1].labels[2],tl.lattice[i+1,j+1].labels[3]], Float64(j)+0.4, Float64(i)+0.5))
            end
        end
    end
    
    for j in 1:n
        if j > 1
            push!(tensors,PlanarTensor(t10,[tl.lattice[n,j].labels[2]], Float64(j)-0.3, Float64(n+1)))
        end
        push!(tensors,PlanarTensor(t01,[tl.lattice[n,j].labels[5]], Float64(j), Float64(n+1)))
        push!(tensors,PlanarTensor(t11,[tl.lattice[n,j].labels[6]], Float64(j)+ 0.3, Float64(n+1)))
    end

    for i in 1:n
        if i > 1
            push!(tensors,PlanarTensor(t11,[tl.lattice[i,n].labels[7]], Float64(n+1), Float64(i)-0.3))
        end
        push!(tensors,PlanarTensor(t01,[tl.lattice[i,n].labels[8]], Float64(n+1), Float64(i)))
        if i < n
            push!(tensors,PlanarTensor(t11,[tl.lattice[i,n].labels[6]], Float64(n+1), Float64(i)+0.3))
        end
    end

    return PlanarTensorNetwork(tensors,max_label)
end