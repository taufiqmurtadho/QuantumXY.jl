using LinearAlgebra

#Constructing the Hamiltonian by its action to Pauli-Z basis states
function method1_hamiltonian(h::Float64, J::Float64, N::Int64)
    #Generating all possible bit sequence of length N
    mat = zeros(2^N, 2^N)
    bitsequence = reverse([reverse(digits(UInt(i), base =2, pad = N)) for i in 0:2^N-1])
    for sequence in bitsequence
        if sequence[N] != sequence[1]
            newsequence = copy(sequence)
            newsequence[N], newsequence[1] = sequence[1], sequence[N]
            mat[findfirst(isequal(sequence), bitsequence), findfirst(isequal(newsequence), bitsequence)]+= -J
        end
        for i in 1:N-1
            if sequence[i]!=sequence[i+1]
                newsequence = copy(sequence)
                newsequence[i], newsequence[i+1] = sequence[i+1], sequence[i]
                mat[findfirst(isequal(sequence), bitsequence), findfirst(isequal(newsequence), bitsequence)]+= -J
            end
        end
    end
    for i in 1:2^N
        mat[i,i] = -h*sum(bitsequence[i]) + h*(N-sum(bitsequence[i]))
    end
    return mat
end

#Constructing the Hamiltonian by kronecker product
function method2_hamiltonian(h::Float64, J::Float64, N::Int64)
    pauliZ = [1. 0; 0 -1.]
    pauliPlus = [0 0;  1. 0]
    pauliMinus = [0 1.; 0 0]
    mat = zeros(2^N, 2^N)
    for i in 1:N
        if i==1
            mat+=-h*kron(pauliZ,Matrix(1.0I,2^(N-1),2^(N-1)))
            hopping = -J*kron(kron(pauliPlus, pauliMinus), Matrix(1.0I, 2^(N-2), 2^(N-2)))
            mat+= hopping+transpose(hopping)
        elseif i == N
            mat+= -h*kron(Matrix(1.0I, 2^(N-1), 2^(N-1)), pauliZ)
            hopping = -J*kron(kron(pauliMinus, Matrix(1.0I, 2^(N-2), 2^(N-2))),pauliPlus)
            mat+= hopping+transpose(hopping)
        else
            mat+=-h*kron(kron(Matrix(1.0I, 2^(i-1), 2^(i-1)), pauliZ), Matrix(1.0I, 2^(N-i), 2^(N-i)))
            hopping = -J*kron(kron(Matrix(1.0I, 2^(i-1), 2^(i-1)), kron(pauliPlus, pauliMinus)), Matrix(1.0I, 2^(N-i-1), 2^(N-i-1)))
            mat+= hopping+transpose(hopping)
        end
    end
    return mat
end


