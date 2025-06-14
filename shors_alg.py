# filename: shors_alg.py
# authors: Amanda Curtis and Arthur Diep-Nguyen

from qiskit.circuit import QuantumCircuit, QuantumRegister, AncillaRegister
from qiskit.quantum_info import Statevector
from qiskit.circuit.library import QFT
#import matplotlib.pyplot as plt

import numpy as np
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
            CU_a |x>_1 |y>_n
        that returns |x>_1 |a*y mod N>_n when x=1,
        where n = 2*ceil(log2(N))
    """
    #TODO see Figure 7
    #TODO don't forget to use to_gate() when composing subcircuits
    pass

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
    quantum_register = QuantumRegister(size=n, name ='x')
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
    quantum_register = QuantumRegister(size=n, name ='x')
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
    phi_b_register = QuantumRegister(size=n, name='phi(b)')

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
    phi_b_register = QuantumRegister(size=n, name='phi(b)')

    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Setting up the circuit 
    c_phi_sub_a = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate((phi_b_register)):
        c_phi_sub_a.cp(np.pi * -a / (2**(n-1-idx)), zero_register[0], q)

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
    phi_b_register = QuantumRegister(size=n, name='phi(b)')
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
    phi_b_register = QuantumRegister(size=n, name='phi(b)')
    # This register is for the extra qubit required
    # Naming to avoid confusion with the function input, a
    zero_register = AncillaRegister(size=1, name="zero")

    # Setting up the circuit 
    cc_phi_sub_a = QuantumCircuit(control_register, phi_b_register, zero_register, name='cc_phi_add_a')

    # Building cc_P_n(a) by making a phase gate p
    # for each qubit 
    for idx, q in enumerate((phi_b_register)):
        cc_phi_sub_a.cp(np.pi * -a / (2**(n-1-idx)), control_register[:], q)

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
    phi_b_register = QuantumRegister(size=n, name='phi(b)')
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
    adder_mod_N.compose(QFT(n).inverse(), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 4/13 controlled NOT (.cx)
    adder_mod_N.cx(phi_b_register[n-1], zero_register[0])

    adder_mod_N.barrier()

    # gate 5/13 qft
    adder_mod_N.compose(QFT(n), phi_b_register[:], inplace=True)

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
    adder_mod_N.compose(QFT(n).inverse(), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 9/13 NOT (aka x)
    adder_mod_N.x(phi_b_register[n-1])

    adder_mod_N.barrier()

    # gate 10/13 controlled NOT (cx)
    adder_mod_N.cx(phi_b_register[n-1], zero_register[0])

    adder_mod_N.barrier()

    # gate 11/13 NOT (aka X)
    adder_mod_N.x(phi_b_register[n-1])

    adder_mod_N.barrier()

    # gate 12/13 qft 
    adder_mod_N.compose(QFT(n), phi_b_register[:], inplace=True)

    adder_mod_N.barrier()

    # gate 13/13 doubly controlled adder 
    cc_sub_a = cc_subtractor(a,N)
    adder_mod_N.compose(cc_sub_a, qubits=[*control_register[:], *phi_b_register[:], *zero_register[:]], inplace=True

    return adder_mod_N

##########################
    

def c_mult_mod(a: int, N: int) -> QuantumCircuit:
    #TODO controlled, see Section 2.3 and Figure 6
    pass

#TODO we also need to write a controlled swap circuit, since we are
# probably not allowed to use the built in swap gate :(

########################################################################
# Other helper functions
########################################################################

#TODO anything that is not a circuit can go here
