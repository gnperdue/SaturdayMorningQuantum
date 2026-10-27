"""The dial demo's code, written the compact way with NumPy arrays.

This file does the same calculations as dial_demo.py, which the lecture imports. Instead
of scoring one gate sequence at a time, it hands NumPy a stack of 24 matrices (one for
each Clifford gate) and scores them all in one operation. That is a little faster, and
much harder to read: start with dial_demo.py, then compare the two.
"""

import numpy as np
from qiskit.synthesis import gridsynth_rz

GATE_MATRICES = {
    'h': np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2),
    's': np.diag([1, 1j]),
    'sdg': np.diag([1, -1j]),
    't': np.diag([1, np.exp(1j * np.pi / 4)]),
    'tdg': np.diag([1, np.exp(-1j * np.pi / 4)]),
    'x': np.array([[0, 1], [1, 0]], dtype=complex),
    'y': np.array([[0, -1j], [1j, 0]]),
    'z': np.diag([1, -1]).astype(complex),
}


def gate_product(sequence):
    """The 2x2 matrix of a list of gate names, applied in order (the first gate acts first)."""
    matrix = np.eye(2, dtype=complex)
    for name in sequence:
        if name not in GATE_MATRICES:
            raise ValueError(f"no matrix for the gate {name!r}")
        matrix = GATE_MATRICES[name] @ matrix
    return matrix


def t_count(sequence):
    """How many T gates (T or its inverse, tdg) a list of gates uses."""
    return sum(name in ('t', 'tdg') for name in sequence)


def sequence_as_qasm(sequence, qubit='q[0]'):
    """A list of gates as lines of QASM."""
    return '\n'.join(f'{name} {qubit};' for name in sequence)


def clifford_sequences():
    """The 24 different one-qubit Clifford gates, each as a short list of H and S gates."""
    found = {}
    frontier = [[]]
    while frontier:
        next_frontier = []
        for sequence in frontier:
            matrix = gate_product(sequence)
            pivot = matrix.flat[np.argmax(np.abs(matrix.flat) > 1e-9)]
            # ignore the overall phase
            key = tuple(np.round((matrix * abs(pivot) / pivot).flatten(), 6))
            if key not in found:
                found[key] = sequence
                next_frontier += [sequence + ['h'], sequence + ['s']]
        frontier = next_frontier
    return list(found.values())


def best_by_t_budget(score, max_t):
    """For each T budget 0..max_t, the lowest-scoring (score, gate list) within that budget.

    Every one-qubit operator made of H, S and T can be written in exactly one way as
    a Clifford gate, then blocks of 'T H' or 'T H S', then possibly one more T, and that
    form uses the fewest T gates possible (Matsumoto and Amano, arXiv:0806.3834). So
    trying every such form tries every operator once. `score` takes an array of 2x2
    matrices and returns one number for each (lower is better).
    """
    cliffords = clifford_sequences()
    clifford_matrices = np.array([gate_product(c) for c in cliffords])
    best = {}
    layer = [([], np.eye(2, dtype=complex))]   # the gates after the Clifford, and their matrix
    for k in range(max_t + 1):
        for middle, matrix in layer:
            for last, n_t in (([], k), (['t'], k + 1)):
                if n_t > max_t:
                    continue
                scores = score(gate_product(last) @ matrix @ clifford_matrices)
                i = int(np.argmin(scores))
                if n_t not in best or scores[i] < best[n_t][0]:
                    best[n_t] = (scores[i], cliffords[i] + middle + last)
        if k < max_t:
            layer = [(middle + block, gate_product(block) @ matrix)
                     for middle, matrix in layer for block in (['t', 'h'], ['t', 'h', 's'])]
    return best_within_budget(best, max_t)


def best_within_budget(best, max_t):
    """From the best (score, gate list) for each exact T count, the best within each budget."""
    within_budget, record = [], None
    for n_t in range(max_t + 1):
        if record is None or best[n_t][0] < record[0]:
            record = best[n_t]
        within_budget.append(record)
    return within_budget


def chance_of_one(matrices):
    """The chance of measuring 1 after H, the gates, H, starting from |0>."""
    hadamard = GATE_MATRICES['h']
    return np.abs((hadamard @ matrices @ hadamard)[..., 1, 0]) ** 2


def best_coins(target, max_t):
    """For each T budget, (chance of 1, gate list) for the coin closest to the target chance."""
    found = best_by_t_budget(lambda m: np.abs(chance_of_one(m) - target), max_t)
    return [(float(chance_of_one(gate_product(sequence))), sequence) for _, sequence in found]


def rotation_error_degrees(matrices, angle):
    """How far gates are from a turn of `angle` radians about z, as the leftover turn in degrees."""
    target = np.diag([np.exp(-0.5j * angle), np.exp(0.5j * angle)])
    # what remains after undoing the target, with the overall phase dropped
    leftover = target.conj().T @ matrices
    leftover = leftover / np.sqrt(np.linalg.det(leftover))[..., None, None]
    # leftover = cos(a) I - i sin(a) (a Pauli direction): a turn by 2a; this form stays accurate
    # for tiny turns, where 2 arccos(|trace| / 2) would be swamped by rounding
    sin_a = np.sqrt(np.abs(leftover[..., 0, 1]) ** 2 + leftover[..., 0, 0].imag ** 2)
    return np.degrees(2 * np.arcsin(np.clip(sin_a, 0, 1)))


def best_rotations(angle, max_t):
    """For each T budget, the gate list closest to a turn of `angle` radians about z."""
    found = best_by_t_budget(lambda m: rotation_error_degrees(m, angle), max_t)
    return [sequence for _, sequence in found]


def qiskit_rotation(angle, epsilon):
    """The gate list Qiskit's compiler (Ross-Selinger) gives for a turn of `angle` about z.

    Qiskit 2.3.1 crashes on many angles when epsilon is larger than about 1e-4 (a bug in
    its gridsynth_rz). The angles in the lecture are tested; if you try others, keep
    epsilon at 1e-6 or smaller.
    """
    circuit = gridsynth_rz(angle, epsilon)
    return [instruction.operation.name for instruction in circuit.data]


def worst_drift(step, step_angle):
    """The largest error in the chance of Left over 20 steps, computed exactly (no shots)."""
    hadamard = GATE_MATRICES['h']
    step_matrix = gate_product(step)
    return max(abs(abs((hadamard @ np.linalg.matrix_power(step_matrix, k) @ hadamard)[1, 0]) ** 2
                   - np.sin(k * step_angle / 2) ** 2)
               for k in range(21))
