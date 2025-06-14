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

# Adder and Subtractor
##############
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
    quantum_register = QuantumRegister(size=n+1, name ='x')
    phi_add_a = QuantumCircuit(quantum_register, name="phi_add_a")

    # Building P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate(reversed(quantum_register)):
        phi_add_a.p(np.pi * a / (1 << idx), q)

    return phi_add_a

def subtractor(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.2, used in figure 5's construction

    # "Solving" for n, necessary to determine number of qubits 
    n = math.ceil(math.log2(N))

    # Setting up Quantum Register 
    quantum_register = QuantumRegister(size=n+1, name ='x')
    phi_sub_a = QuantumCircuit(quantum_register, name="phi_sub_a")

    # Building rev(P_n(a)) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate(reversed(quantum_register)):
        phi_sub_a.p(np.pi * -a / (1 << idx), q)

    return phi_sub_a

# Controlled adder and subtractor 
# Needed to Build Modular Adder Gate
#################

# c adder gate
def c_adder(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits 
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')

    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Setting up the circuit 
    c_phi_add_a = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate(reversed(phi_b_register)):
        c_phi_add_a.cp(np.pi * a / (1 << idx), zero_register[0], q)

    return c_phi_add_a

# c subtractor gate
def c_subtractor(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits 
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')

    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Setting up the circuit 
    c_phi_sub_a = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate((phi_b_register)):
        c_phi_sub_a.cp(np.pi * -a / (2**(n-idx)), zero_register[0], q)

    return c_phi_sub_a

# Doubly-controlled adder and subtractor 
# Needed to Build Modular Adder Gate
#################

# cc adder gate
def cc_adder(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits 
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')
    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Setting up the circuit 
    cc_phi_add_a = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate(reversed(phi_b_register)):
        cc_phi_add_a.mcp(np.pi * a / (1 << idx), control_register[:], q)

    return cc_phi_add_a

# cc subtractor gate 
def cc_subtractor(a: int, N: int) -> QuantumCircuit:
    #TODO see Section 2.1 and Figure 3
    # Must first "solve" for n, set up registers needed

    # "Solving" for n, necessary to determine number of qubits 
    n = math.ceil(math.log2(N))

    # Setting up Quantum Registers
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')
    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Setting up the circuit 
    cc_phi_sub_a = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate((phi_b_register)):
        cc_phi_sub_a.cp(np.pi * -a / (2**(n-idx)), control_register[:], q)

    return cc_phi_sub_a



def cc_adder_mod(a: int, N: int) -> QuantumCircuit:
    # Copying outline of code from our .py file 

def adder_mod(a: int, N: int) -> QuantumCircuit:
    #TODO doubly controlled, see Section 2.2 and Figure 5
   
   # "Solving" for n, necessary to determine number of qubits 
    n = math.ceil(math.log2(N))

    # Makking the control register, which requires two qubits
    control_register = QuantumRegister(size=2, name='c')
    # This is the register where most of the work happens
    phi_b_register = QuantumRegister(size=n+1, name='phi(b)')
    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Creating the actual circuit 
    # naming it 'adder_mod_N'
    adder_mod_N = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_adder_mod_N')

    # gate 1/13 - cc_adder with a
    # ensuring cleaner code by naming gate and specifying a, N
    cc_adder_aN = cc_adder(a,N)
    adder_mod_N.compose(cc_adder_aN, qubits=[*control_register[:], *phi_b_register[:], *zero_register[:]], inplace=True)

    adder_mod_N.barrier()

    # gate 2/13 inverse adder (aka subtractor) with a = N
    sub_N = subtractor(N,N)
    adder_mod_N.compose(sub_N, qubits=[*phi_b_register[:]], inplace=True)

    adder_mod_N.barrier()

    # gate 3/13 qft_inv
    adder_mod_N.compose(QFT(n+1).inverse(), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 4/13 controlled NOT (.cx)
    adder_mod_N.cx(phi_b_register[n], zero_register[0])

    adder_mod_N.barrier()

    # gate 5/13 qft
    adder_mod_N.compose(QFT(n+1), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 6/13 controlled adder with a = N 
    c_add_N = c_adder(N,N)
    adder_mod_N.compose(c_add_N, qubits=[*control_register[:], *phi_b_register[:], *zero_register[:]], inplace=True)

    adder_mod_N.barrier()

    # gate 7/13 doubly controlled adder inverse (aka doubly controlled subtractor) for a
    cc_sub_a = cc_subtractor(a,N)
    adder_mod_N.compose(cc_sub_a, qubits=[*control_register[:], *phi_b_register[:], *zero_register[:]], inplace=True)

    adder_mod_N.barrier()

    # gate 8/13 qft inverse on phi_b
    adder_mod_N.compose(QFT(n+1).inverse(), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 9/13 NOT (aka x)
    adder_mod_N.x(phi_b_register[n])

    adder_mod_N.barrier()

    # gate 10/13 controlled NOT (cx)
    adder_mod_N.cx(phi_b_register[n], zero_register[0])

    adder_mod_N.barrier()

    # gate 11/13 NOT (aka X)
    adder_mod_N.x(phi_b_register[n])

    adder_mod_N.barrier()

    # gate 12/13 qft 
    adder_mod_N.compose(QFT(n+1), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 13/13 doubly controlled adder 
    cc_sub_a = cc_subtractor(a,N)
    adder_mod_N.compose(cc_sub_a, qubits=[*control_register[:], *phi_b_register[:], *zero_register[:]], inplace=True

    return adder_mod_N

##########################
    

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
