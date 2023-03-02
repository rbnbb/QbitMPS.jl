using FFTW

"Compute Quantum Fourier Transfor of number state using FFT"
function fft_statevector(number::Int, numbits::Int)
    @assert number < 2^numbits "Cannot represent $number in $numbits bits"
    input_state = zeros(2^numbits)
    # we count from the end since the QFT has a positive phase,
    # but FFTW::fft() implements DFT with negative phase
    input_state[2^numbits-number+1] = 1
    return fft(input_state)/sqrt(2^numbits)
end

"Compute Quantum Fourier Transfor of number state using MPS"
function mps_statevector(number::Int, numbits::Int; max_bond_dimension=128)
    @assert number < 2^numbits "Cannot represent $number in $numbits bits"
    input_state = [string(x) for x in QbitMPS.int2binary(number, numbits)]
    qft_circuit = QbitMPS.quantum_fourier_circuit(numbits)
    psi = simulate_circuit(numbits, qft_circuit; initial_state=input_state, maxdim=max_bond_dimension)
    return QbitMPS.mps2statevector(psi)
end
