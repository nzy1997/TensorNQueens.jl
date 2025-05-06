struct PlanarTensor{T <: Integer, T2}
    tensor::Array{T2}
    labels::Vector{T}
    x::Float64
    y::Float64
end
struct PlanarTensorNetwork{T}
    tensors::Vector{PlanarTensor{Int, T}}
    max_label::Int
end