# filename: shors_alg.py
# authors: Amanda Curtis and Arthur Diep-Nguyen

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.quantum_info import Statevector
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
        Positive integer
    N : int
        Positive integer strictly greater than `a`

    Returns
    -------
    QuantumCircuit
        Controlled circuit CU_a given by
            CU_a |c>_1 |x>_n
        that returns |c>_1 |a*x mod N>_n when c=1,
        where n = ceil(log2(N))
    """
    #TODO see Figure 7
    #TODO don't forget to use to_gate() when composing subcircuits

    # Number of bits required to represent N
    n = math.ceil(math.log2(N))

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
    qc.compose(c_mult_mod(a,N),
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

    qc.compose(c_mult_mod(a,N).inverse(),
                qubits=[control_qr, input_qr, ancilla],
                inplace=True)

    return qc

########################################################################
# Subcircuits
########################################################################

def adder(a: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    pass

def cc_adder_mod(a: int, N: int) -> QuantumCircuit:
    #TODO doubly controlled, see Section 2.2 and Figure 5
    pass

def c_mult_mod(a: int, N: int) -> QuantumCircuit:
    #TODO controlled, see Section 2.3 and Figure 6
    pass

########################################################################
# Other helper functions
########################################################################

#TODO anything that is not a circuit can go here
