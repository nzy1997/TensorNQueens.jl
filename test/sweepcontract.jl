using TensorNQueens
using Test
using TensorNQueens:generate_planar_tensor_network
using TropicalSweepContractor
using TropicalNumbers
using TropicalSweepContractor.GenericTensorNetworks

@testset "Sweep Contract" begin
    N = 1021
    ft = generate_planar_tensor_network(7,Mod{N,Int})
    @show sweep_contract!(ft,1,1)
end 