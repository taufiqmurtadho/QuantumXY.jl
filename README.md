# QuantumXY.jl
This is a Julia package to simulate a 1D-quantum spin chain on a ring under XY Hamiltonian.

## FUNCTION
```method1_hamiltonian(h, J, N)```: Constructing the Hamiltonian by its action toward Pauli-Z basis states
```method2_hamiltonian(h, J, N)```: Constructing the Hamiltonian by direct kronecker product

## INPUT
h::Float64 (magnetic field strength), J::Float64 (nearest-neighbor interaction strength), N::Int64(system's size)

## OUTPUT 
$2^N$ x $2^N$ array representing Hamiltonian matrix in Pauli-Z basis
