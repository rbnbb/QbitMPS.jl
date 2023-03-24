using StateVecSim

function check_statevector(circuit, numqubits::Integer)
    statevector = zero_qreg(numqubits)
    StateVecSim.apply!(circuit, statevector)
    return statevector.real + 1im * statevector.imag
end
