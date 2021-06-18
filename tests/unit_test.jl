using Test
using LinearAlgebra
using QuantumXY

@testset "Check methods equivalence" begin
	for n in 2:9
		@test method1_hamiltonian(0.5,1.0,n) == method2_hamiltonian(0.5,1.0,n)
	end
end

@testset "Check ground state J =0" begin
	#h > 0 
	eig1 = eigen(method1_hamiltonian(0.5,0.0,3))
	gs_index = findfirst(isequal(min(eig1.values ...)), eig1.values)
	@test eig1.vectors[:,gs_index] == [1., 0., 0., 0., 0., 0., 0., 0.]
	# h < 0
	eig2 = eigen(method2_hamiltonian(-0.5,0.0,3))
	gs_index = findfirst(isequal(min(eig2.values ...)), eig2.values)
	@test eig2.vectors[:,gs_index] == [0., 0., 0., 0., 0., 0., 0., 1.]
end
