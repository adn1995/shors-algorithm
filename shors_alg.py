# filename: shors_alg.py
# authors: Amanda Curtis and Arthur Diep-Nguyen

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.quantum_info import Statevector
from qiskit.circuit.library import QFT
#import matplotlib.pyplot as plt

import math

########################################################################
# Main oracle function
########################################################################

def oracle(a: int, N: int) -> QuantumCircuit:
    """Returns the oracle from Shor's algorithm.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator U given by
            U |c>_1 |x>_n =
                |c>_1 |a*x mod N>_n if c==1 and x<N
                |c>_1 |x>_n         otherwise
        where n = ceil(log2(N))
    """
    # See Figure 7 of N&C

    # Number of bits required to represent N
    n = math.ceil(math.log2(N))

    # TODO/WARNING: THIS FUNCTION SHOULD CHECK IF x < N

    # One control qubit
    control_qr = QuantumRegister(1, name="c")

    # n qubits for |x>, the n-qubit input for CMULT(a)MOD(N)
    input_qr = QuantumRegister(n, name="x")

    # n+2 qubits
    # n+1 qubits are for (n+1)-qubit input |b> for CMULT(a)MOD(N)
    # The last qubit is the ancillary qubit
    ancilla = AncillaRegister(n+2, name="a")

    qc = QuantumCircuit(control_qr, input_qr, ancilla, name="oracle")

    # Controlled multiplier
    qc.compose(c_mult_mod(a,N).to_gate(),
                qubits=[control_qr, input_qr, ancilla],
                inplace=True)


    # CMULT(a)MOD(N) has the following output wires
    #   |c>             1 control qubit
    #   |x>             n qubits
    #   |b+a*x mod N>   n+1 qubits
    #   |0>             1 ancillary qubit from ADD(a)MOD(N)
    #
    # We want to controlled-swap |x> with |b+a*x mod N>,
    # controlled by |c>.
    #
    # Since b=0, we have |b+a*x mod N> = |ax mod N>, which requires
    # only n qubits, not n+1 qubits, so we can ignore the most
    # significant qubit of |b+a*x mod N>.
    #
    # In other words, we want to controlled-swap the n qubits of |x>
    # with the n least significant qubits of |b+a*x mod N>.
    for i in range(n):
        qc.cswap(control_qubit=control_qr,
            target_qubit1=input_qr[i],
            target_qubit2=ancilla[i+1])

    # Inverse controlled multiplier

    qc.compose(c_mult_mod(a,N).inverse().to_gate(),
                qubits=[control_qr, input_qr, ancilla],
                inplace=True)

    return qc

########################################################################
# Subcircuits
########################################################################

def adder(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    # Must first "solve" for n, set up registers needed
    # Based off of cited paper, QC Bootcamp Problem Session 2's implementation of
    # Draper's adder circuit, and class lecture on 22 May 2025

    # We utilize the corollary that relates the QFT, A_k (Draper),
    # and this (controlled) phase gate P_n(a) (phiADD(a))

    # "Solving" for n
    n = math.ceil(math.log2(N))

    # Setting up Quantum Register
    quantum_register = QuantumRegister(size=n, name ='x')
    phi_add_a = QuantumCircuit(quantum_register, name="phi_add_a")

    # Building P_n(a) by making a phase gate p
    # for each qubit
    for idx, q in enumerate(reversed(quantum_register)):
        phi_add_a.p(math.pi * a / (1 << idx), q)

    return phi_add_a


def cc_adder_mod(a: int, N: int) -> QuantumCircuit:
    #TODO doubly controlled, see Section 2.2 and Figure 5
    pass

def c_mult_mod(a: int, N: int) -> QuantumCircuit:
    """Returns the CMULT(a)MOD(N) circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Controlled multiplier circuit that takes 3 inputs
            |c>_1     control qubit
            |x>_n
            |b>_n
        with the output
            |c> |x> |b>             if c==0
            |c> |x> |b+a*x mod N>   if c==1
    """
    # See Section 2.3 and Figure 6 of N&C

    # Number of bits required to represent N
    n = math.ceil(math.log2(N))

    # One control qubit
    control_qr = QuantumRegister(1, name="c")

    # n qubits for |x>, the n-qubit input for CMULT(a)MOD(N)
    input1_qr = QuantumRegister(n, name="x")

    # n+1 qubits for |b>
    # Even though b is an n-bit number, we need n+1 qubits
    # to account for overflow
    input2_qr = QuantumRegister(n+1, name="b")

    # 1 ancilla qubit used by modular adder
    ancilla = AncillaRegister(1, name="a")

    qc = QuantumCircuit(control_qr, input1_qr, input2_qr, ancilla,
                        name="CMULT({})MOD({})".format(str(a),str(N)))

    # QFT circuit for n+1 qubit register
    qft = QFT(n+1)

    # Apply QFT to |b> register
    qc.compose(qft.to_gate(), input2_qr, inplace=True)

    # n doubly-controlled modular adders
    for i in range(n):
        qc.compose(cc_adder_mod((2**i)*a, N).to_gate(),
                    qubits=[control_qr,
                        input1_qr[i],
                        *input2_qr,
                        ancilla],
                    inplace=True)

    # Apply inverse QFT to |b> register
    qc.compose(qft.inverse().to_gate(), input2_qr, inplace=True)

    return qc

########################################################################
# Other helper functions
########################################################################

#TODO anything that is not a circuit can go here
