using TensorNQueens
using Test
using TensorNQueens:generate_planar_tensor_network
using TropicalSweepContractor
using TropicalNumbers

@testset "Sweep Contract" begin
    ft = generate_planar_tensor_network(5,TropicalAndOr)
    @show sweep_contract!(ft,10,10)
end