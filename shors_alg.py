# filename: shors_alg.py
# authors: Amanda Curtis and Arthur Diep-Nguyen
"""This module contains functions that construct quantum circuits used
in Shor's algorithm.

The implementations of these circuits come from the paper "Circuit for
Shor's Algorithm using 2n+3 qubits" by Stephane Beauregard.
"""

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.circuit.library import QFT
#import matplotlib.pyplot as plt

import math

########################################################################
########################################################################
# Main oracle function
########################################################################
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
    # See Figure 7 of "Circuit for Shor's Algorithm using 2n+3 qubits"

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
    qc.compose(c_mult_mod_inv(a,N).to_gate(),
                qubits=[control_qr, input_qr, ancilla],
                inplace=True)

    return qc

########################################################################
########################################################################
# Subcircuits
########################################################################
########################################################################

########################################################################
# Adder and subtractor
########################################################################

def adder(a: int, N: int) -> QuantumCircuit:
    """Returns the uncontrolled adder circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes in the state
            |phi(b)>_(n+1)
        and outputs the state
            |phi(a+b)>_(n+1)
    """
    # See Section 2.1 and Figure 3 of "Circuit for Shor's Algorithm
    # using 2n+3 qubits"

    # Must first "solve" for n, set up registers needed
    # Based off of cited paper, QC Bootcamp Problem Session 2's
    # implementation of Draper's adder circuit, and class lecture on
    # 22 May 2025

    # We utilize the corollary that relates the QFT, A_k (Draper),
    # and this (controlled) phase gate P_n(a) (phiADD(a))

    # "Solving" for n
    n = math.ceil(math.log2(N))

    # Setting up Quantum Register
    quantum_register = QuantumRegister(size=n+1, name ='x')
    phi_add_a = QuantumCircuit(quantum_register, name="phi_add_a")

    # Building P_n(a) by making a phase gate p
    # for each qubit
    for idx, q in enumerate(reversed(quantum_register)):
        phi_add_a.p(math.pi * a / (1 << idx), q)

    return phi_add_a

def subtractor(a: int, N: int) -> QuantumCircuit:
    """Returns the uncontrolled subtractor circuit, i.e. the inverse
    of the adder circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes in the state
            |phi(b)>_(n+1)
        and outputs the state
            |phi(b-a)>_(n+1)
    """
    # See Section 2.2 of "Circuit for Shor's Algorithm using 2n+3
    # qubits"
    # Used in the construction of the modular adder from Figure 5

    # "Solving" for n, necessary to determine number of qubits
    n = math.ceil(math.log2(N))

    # Setting up Quantum Register
    quantum_register = QuantumRegister(size=n+1, name ='x')
    phi_sub_a = QuantumCircuit(quantum_register, name="phi_sub_a")

    # Building rev(P_n(a)) by making a phase gate p
    # for each qubit
    for idx, q in enumerate(reversed(quantum_register)):
        phi_sub_a.p(math.pi * -a / (1 << idx), q)

    return phi_sub_a

########################################################################
# Singly-controlled adder and subtractor
# Needed to build modular adder gate
########################################################################

def c_adder(a: int, N: int) -> QuantumCircuit:
    """Returns the singly-controlled adder circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes the input state
            |c>_1 |phi(b)>_(n+1)
        and outputs the state
            |c>_1 |phi(a+b)>_(n+1)  if c==1
            |c>_1 |phi(b)>_(n+1)    if c==0
    """
    # See Section 2.1 and Figure 3 of "Circuit for Shor's Algorithm
    # using 2n+3 qubits"
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=1, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')

    # Setting up the circuit
    c_phi_add_a = QuantumCircuit(control_register,
                                    phi_b_register,
                                    name='c_phi_add_a')

    # Building c_P_n(a) by making a phase gate p
    # for each qubit
    for idx, q in enumerate(reversed(phi_b_register)):
        c_phi_add_a.cp(math.pi * a / (1 << idx), control_register, q)

    return c_phi_add_a

def c_subtractor(a: int, N: int) -> QuantumCircuit:
    """Returns the singly-controlled subtractor circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes the input state
            |c>_1 |phi(b)>_(n+1)
        and outputs the state
            |c>_1 |phi(b-a)>_(n+1)  if c==1
            |c>_1 |phi(b)>_(n+1)    if c==0
    """
    # See Section 2.1 and Figure 3 of "Circuit for Shor's Algorithm
    # using 2n+3 qubits"
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=1, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')

    # Setting up the circuit
    c_phi_sub_a = QuantumCircuit(control_register,
                                    phi_b_register,
                                    name='c_phi_sub_a')

    # Building c_P_n(a) by making a phase gate p
    # for each qubit
    for idx, q in enumerate((phi_b_register)):
        c_phi_sub_a.cp(math.pi * -a / (2**(n-idx)), control_register, q)

    return c_phi_sub_a

########################################################################
# Doubly-controlled adder and subtractor
# Needed to build modular adder gate
########################################################################

def cc_adder(a: int, N: int) -> QuantumCircuit:
    """Returns the doubly-controlled adder circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes the input state
            |c_1 c_2>_2 |phi(b)>_(n+1)
        and outputs the state
            |c_1 c_2>_2 |phi(a+b)>_(n+1)  if c_1==1 and c_2==1
            |c_1 c_2>_2 |phi(b)>_(n+1)    if otherwise
    """
    # See Section 2.1 and Figure 3 of "Circuit for Shor's Algorithm
    # using 2n+3 qubits"
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')

    # Setting up the circuit
    cc_phi_add_a = QuantumCircuit(control_register,
                                    phi_b_register,
                                    name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit
    for idx, q in enumerate(reversed(phi_b_register)):
        cc_phi_add_a.mcp(math.pi * a / (1 << idx), control_register[:], q)

    return cc_phi_add_a

def cc_subtractor(a: int, N: int) -> QuantumCircuit:
    """Returns the doubly-controlled subtractor circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes the input state
            |c_1 c_2>_2 |phi(b)>_(n+1)
        and outputs the state
            |c_1 c_2>_2 |phi(b-a)>_(n+1)  if c_1==1 and c_2==1
            |c_1 c_2>_2 |phi(b)>_(n+1)    if otherwise
    """
    # See Section 2.1 and Figure 3 of "Circuit for Shor's Algorithm
    # using 2n+3 qubits"
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')

    # Setting up the circuit
    cc_phi_sub_a = QuantumCircuit(control_register,
                                    phi_b_register,
                                    name='cc_phi_sub_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit
    for idx, q in enumerate((phi_b_register)):
        cc_phi_sub_a.mcp(math.pi * -a / (2**(n-idx)), control_register[:], q)

    return cc_phi_sub_a

########################################################################
# Doubly-controlled modular adder (and its inverse)
########################################################################

def cc_adder_mod(a: int, N: int) -> QuantumCircuit:
    """Returns the doubly-controlled modular adder circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes the input state
            |c_1 c_2>_2 |phi(b)>_(n+1)
        and outputs the state
            |c_1 c_2>_2 |phi(a+b) mod N>_(n+1)  if c_1==1 and c_2==1
            |c_1 c_2>_2 |phi(b)>_(n+1)          otherwise
    """
    # See Section 2.2 and Figure 5 of "Circuit for Shor's Algorithm
    # using 2n+3 qubits"

    # "Solving" for n, necessary to determine number of qubits
    n = math.ceil(math.log2(N))

    # Making the control register, which requires two qubits
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')
    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Creating the actual circuit
    # naming it 'adder_mod_N'
    adder_mod_N = QuantumCircuit(control_register,
                                    phi_b_register,
                                    zero_register,
                                    name='cc_adder_mod_N')

    # Subcircuits that we will use later
    c_add_N = c_adder(N,N).to_gate()
    cc_adder_aN = cc_adder(a,N).to_gate()
    sub_N = subtractor(N,N).to_gate()
    cc_sub_aN = cc_subtractor(a,N).to_gate()

    # gate 1/13 - cc_adder with a
    adder_mod_N.compose(cc_adder_aN,
                        qubits=[*control_register[:],
                                *phi_b_register[:]],
                        inplace=True)

    #adder_mod_N.barrier()

    # gate 2/13 inverse adder (aka subtractor) with a = N
    adder_mod_N.compose(sub_N, qubits=[*phi_b_register[:]], inplace=True)

    #adder_mod_N.barrier()

    # gate 3/13 qft_inv
    adder_mod_N.compose(QFT(n+1).inverse(), phi_b_register[:], inplace=True)

    #adder_mod_N.barrier()

    # gate 4/13 controlled NOT (.cx)
    adder_mod_N.cx(phi_b_register[n], zero_register[0])

    #adder_mod_N.barrier()

    # gate 5/13 qft
    adder_mod_N.compose(QFT(n+1), phi_b_register[:], inplace=True)

    #adder_mod_N.barrier()

    # gate 6/13 controlled adder with a = N
    adder_mod_N.compose(c_add_N,
                        qubits=[zero_register[0],
                                *phi_b_register[:]],
                        inplace=True)

    #adder_mod_N.barrier()

    # gate 7/13 doubly controlled adder inverse (aka doubly controlled
    # subtractor) for a
    adder_mod_N.compose(cc_sub_aN,
                        qubits=[*control_register[:],
                                *phi_b_register[:]],
                        inplace=True)

    #adder_mod_N.barrier()

    # gate 8/13 qft inverse on phi_b
    adder_mod_N.compose(QFT(n+1).inverse(), phi_b_register[:], inplace=True)

    #adder_mod_N.barrier()

    # gate 9/13 NOT (aka x)
    adder_mod_N.x(phi_b_register[n])

    #adder_mod_N.barrier()

    # gate 10/13 controlled NOT (cx)
    adder_mod_N.cx(phi_b_register[n], zero_register[0])

    #adder_mod_N.barrier()

    # gate 11/13 NOT (aka X)
    adder_mod_N.x(phi_b_register[n])

    #adder_mod_N.barrier()

    # gate 12/13 qft
    adder_mod_N.compose(QFT(n+1), phi_b_register[:], inplace=True)

    #adder_mod_N.barrier()

    # gate 13/13 doubly controlled adder
    adder_mod_N.compose(cc_sub_aN,
                        qubits=[*control_register[:],
                                *phi_b_register[:]],
                        inplace=True)

    return adder_mod_N

def cc_adder_mod_inv(a: int, N: int) -> QuantumCircuit:
    """Returns the inverse of the doubly controlled modular adder
    circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes the input state
            |c_1 c_2>_2 |phi(b)>_(n+1)
        and outputs the state
            |c_1 c_2>_2 |phi(b-a) mod N>_(n+1)  if c_1==1 and c_2==1
            |c_1 c_2>_2 |phi(b)>_(n+1)          otherwise
    """
    # Number of bits required to represent N
    n = math.ceil(math.log2(N))

    # Two control qubits
    control_qr = QuantumRegister(2, name="c")

    # n+1 qubits for |phi(b)>
    # Even though b is an n-bit number, we need n+1 qubits
    # to account for overflow
    input_qr = QuantumRegister(n+1, name="phi(b)")

    # 1 ancilla qubit
    ancilla = AncillaRegister(1, name="a")

    qc = QuantumCircuit(control_qr, input_qr, ancilla,
                        name="CCphi({})MOD({})inv".format(str(a),str(N)))

    # QFT circuit for n+1 qubit register
    qft = QFT(n+1)

    qc.compose(cc_subtractor(a,N).to_gate(),
                qubits=[*control_qr, *input_qr],
                inplace=True)

    qc.compose(qft.inverse().to_gate(), input_qr, inplace=True)

    qc.x(input_qr[n])

    qc.cx(input_qr[n], ancilla)

    qc.x(input_qr[n])

    qc.compose(qft.to_gate(), input_qr, inplace=True)

    qc.compose(cc_adder(a,N).to_gate(),
                qubits=[*control_qr, *input_qr],
                inplace=True)

    qc.compose(c_subtractor(N,N).to_gate(),
                qubits=[ancilla[0], *input_qr],
                inplace=True)

    qc.compose(qft.inverse().to_gate(), input_qr, inplace=True)

    qc.cx(input_qr[n], ancilla)

    qc.compose(qft.to_gate(), input_qr, inplace=True)

    qc.compose(adder(N,N).to_gate(), input_qr, inplace=True)

    qc.compose(cc_subtractor(a,N).to_gate(),
                qubits=[*control_qr, *input_qr],
                inplace=True)

    return qc

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
        Unitary operator that takes 3 inputs
            |c>_1     control qubit
            |x>_n
            |b>_n
        and outputs the state
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

def c_mult_mod_inv(a: int, N: int) -> QuantumCircuit:
    """Returns the inverse of the CMULT(a)MOD(N) circuit.

    Parameters
    ----------
    a : int
        Positive integer strictly less than `N`
    N : int
        Positive integer, which Shor's algorithm factors

    Returns
    -------
    QuantumCircuit
        Unitary operator that takes 3 inputs
            |c>_1     control qubit
            |x>_n
            |b>_n
        and outputs the state
            |c> |x> |b>             if c==0
            |c> |x> |b-a*x mod N>   if c==1
    """
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

    # n inverses of doubly-controlled modular adder
    for i in range(n):
        qc.compose(cc_adder_mod_inv((2**(n-1-i))*a, N).to_gate(),
                    qubits=[control_qr,
                        input1_qr[n-1-i],
                        *input2_qr,
                        ancilla],
                    inplace=True)

    # Apply inverse QFT to |b> register
    qc.compose(qft.inverse().to_gate(), input2_qr, inplace=True)

    return qc
