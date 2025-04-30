# 3 4 8
# 1 5 9
# 2 6 7
function generate_tensor(T)
    res = ones(T, fill(2,9)...)
    for index in 0:2^9-1
        index = digits(index, base = 2, pad = 9) .+ 1
        tag = false
        if index[5] == 2
            for i in 6:9
                if index[i] == 1
                    tag = true
                end
            end
        end
        for i in 1:4
            if index[i] == 2
                if (index[10-i] == 1) || (index[5] == 2)
                    tag = true
                end
            end
        end
        for i in 6:9
            if index[i] == 2
                if (index[10-i] == 1) && (index[5] == 1)
                    tag = true
                end
            end
        end
        if tag
            res[index...] = zero(T)
        end
    end
    return res
end

function generate01tensors(T)
    return[one(T),zero(T)],[zero(T),one(T)] ,[one(T),one(T)]
end

struct TensorNQ
    labels::Vector{Int}
end

function Base.show(io::IO, t::TensorNQ)
    print(io, "\n$(t.labels[3]) $(t.labels[4]) $(t.labels[8])")
    print(io, "\n$(t.labels[1]) $(t.labels[5]) $(t.labels[9])")
    print(io, "\n$(t.labels[2]) $(t.labels[6]) $(t.labels[7])")
end

# (i-1,j-1) (i-1,j) (i-1,j+1)
# (i,j-1)   (i,j)   (i,j+1)
# (i+1,j-1) (i+1,j) (i+1,j+1)

# 3 4 8
# 1 5 9
# 2 6 7
function get_pos()
    return [(0,-1), (1,-1), (-1,-1), (-1,0),(0,0), (1,0), (1,1),(-1,1),(0,1)]
end

struct TensorNQ_lattice
    lattice::Matrix{TensorNQ}
    pos10::Vector{Int}
    pos01::Vector{Int}
    pos11::Vector{Int}
end

function generate_TensorNQ_lattice(n::Int)
    # res = Matrix{TensorNQ}(undef,n,n)
    res = fill(TensorNQ(fill(-1,9)),n,n)
    label_count = 0

    pos_vec = get_pos()

    pos11 = Int[]
    pos01 = Int[]
    pos10 = Int[]
    for i in 1:n
        for j in 1:n
            index_vec = fill(-1,9)
            for (ind,pos) in enumerate(pos_vec)
                if checkbounds(Bool,res,i+pos[1],j+pos[2])
                    if res[i+pos[1],j+pos[2]].labels[10-ind] != -1
                        index_vec[ind] = res[i+pos[1],j+pos[2]].labels[10-ind]
                    else
                        label_count += 1
                        index_vec[ind] = label_count
                    end
                else
                    label_count += 1
                    index_vec[ind] = label_count
                    if ind < 5
                        push!(pos10,label_count)
                    elseif ind == 9 || ind == 6
                        push!(pos01,label_count)
                    else
                        push!(pos11,label_count)
                    end
                end
            end
            res[i,j] = TensorNQ(index_vec)
        end
    end
    return TensorNQ_lattice(res, pos10,pos01,pos11 ∪ [res[i,j].labels[5] for i in 1:n, j in 1:n])
end

function generate_tensor_network(n::Int,T)
    t9_lattice = generate_TensorNQ_lattice(n)
    lattice,pos10,pos01,pos11 = t9_lattice.lattice,t9_lattice.pos10,t9_lattice.pos01,t9_lattice.pos11
    t9_ixs = getfield.(vec(lattice),:labels)
    t = generate_tensor(T)
    t10,t01,t11 = generate01tensors(T)
    return DynamicEinCode(t9_ixs ∪ [[p] for p in pos10] ∪ [[p] for p in pos01] ∪  [[p] for p in pos11],Int[]),[fill(t,length(t9_ixs))...,fill(t10,length(pos10))...,fill(t01,length(pos01))...,fill(t11,length(pos11))...]
end