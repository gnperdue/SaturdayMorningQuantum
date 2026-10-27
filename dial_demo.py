"""The code behind the lecture's dial demo, written to be read.

The demo builds turns of a qubit out of the few gates a fault-tolerant computer offers
(H, S and T). Every step here is spelled out, so you can follow it with pencil and paper.
Start with `apply_gate` and `chance_of_one`: work out chance_of_one(T) by hand, and you
should get about 0.146, the 14.64% in the lecture's coin table.

How we store things:
- A qubit's state is a list of two amplitudes: [amplitude of |0>, amplitude of |1>].
- A gate is a 2x2 matrix, stored as a list of its two rows: [[a, b], [c, d]].
- A gate sequence is a list of gate names, in the order they act: ['h', 't', 'h'].

Python counts from 0: state[0] is the amplitude of |0>, and gate[1][0] is the entry in
row 1, column 0 (bottom left).

dial_demo_numpy.py does the same calculations the compact way, by having NumPy score 24
matrices at once. It is a little faster and much harder to read: compare the two once
this file makes sense.
"""

import cmath
import math

import matplotlib.pyplot as plt
from qiskit.synthesis import gridsynth_rz

# ----- The gates -----

# Python writes the imaginary unit i as 1j.
ROOT_HALF = 1 / math.sqrt(2)                 # 1/sqrt(2), about 0.707
EIGHTH_TURN = cmath.exp(1j * math.pi / 4)    # e^(i pi/4), the phase T puts on |1>

IDENTITY = [[1, 0], [0, 1]]                  # the "do nothing" gate
H = [[ROOT_HALF, ROOT_HALF], [ROOT_HALF, -ROOT_HALF]]
S = [[1, 0], [0, 1j]]                        # a quarter turn about z
T = [[1, 0], [0, EIGHTH_TURN]]               # an eighth of a turn about z

# Every gate the code may meet, by its QASM name. 'sdg' and 'tdg' undo S and T.
GATE_MATRICES = {
    'h': H,
    's': S,
    'sdg': [[1, 0], [0, -1j]],
    't': T,
    'tdg': [[1, 0], [0, EIGHTH_TURN.conjugate()]],
    'x': [[0, 1], [1, 0]],
    'y': [[0, -1j], [1j, 0]],
    'z': [[1, 0], [0, -1]],
}


# ----- Matrix arithmetic -----

def apply_gate(gate, state):
    """The new state after a gate acts on a state: the matrix times the vector.

    Recall from the lecture: to multiply a matrix by a vector, dot each row of the matrix
    with the vector. For a 2x2 gate g and a state s with two amplitudes:

        | g00  g01 |   | s0 |     | g00 s0 + g01 s1 |
        |          | x |    |  =  |                 |
        | g10  g11 |   | s1 |     | g10 s0 + g11 s1 |

    Here g01 means gate[0][1] (row 0, column 1) and s1 means state[1]. The two lines of
    code below are read straight off the right-hand side.
    """
    new_amplitude_0 = gate[0][0] * state[0] + gate[0][1] * state[1]
    new_amplitude_1 = gate[1][0] * state[0] + gate[1][1] * state[1]
    return [new_amplitude_0, new_amplitude_1]


def multiply(first, second):
    """The product of two 2x2 matrices: the single gate that does `second`, then `first`.

    Each entry of the answer is a row of `first` dotted with a column of `second`:

        | a00  a01 |   | b00  b01 |     | a00 b00 + a01 b10    a00 b01 + a01 b11 |
        |          | x |          |  =  |                                       |
        | a10  a11 |   | b10  b11 |     | a10 b00 + a11 b10    a10 b01 + a11 b11 |

    Mind the order: in a product, the matrix on the right acts first.
    """
    top_left = first[0][0] * second[0][0] + first[0][1] * second[1][0]
    top_right = first[0][0] * second[0][1] + first[0][1] * second[1][1]
    bottom_left = first[1][0] * second[0][0] + first[1][1] * second[1][0]
    bottom_right = first[1][0] * second[0][1] + first[1][1] * second[1][1]
    return [[top_left, top_right], [bottom_left, bottom_right]]


def gate_product(sequence):
    """The single matrix that does a whole list of gates, in order.

    The first gate in the list acts first, so each new gate multiplies on the left:
    for ['h', 't'] the matrix is T x H.
    """
    matrix = IDENTITY
    for name in sequence:
        if name not in GATE_MATRICES:
            raise ValueError(f"no matrix for the gate '{name}'")
        matrix = multiply(GATE_MATRICES[name], matrix)
    return matrix


def t_count(sequence):
    """How many T gates a list of gates uses (T itself, or 'tdg', which undoes a T)."""
    count = 0
    for name in sequence:
        if name == 't' or name == 'tdg':
            count = count + 1
    return count


def sequence_as_qasm(sequence, qubit='q[0]'):
    """A list of gates as QASM, one line per gate: ['h', 't'] becomes 'h q[0];' and 't q[0];'."""
    lines = []
    for name in sequence:
        lines.append(f'{name} {qubit};')
    return '\n'.join(lines)    # one string, with a line break between the lines


# ----- The warm-up: a coin -----

def chance_of_one(gates_matrix):
    """The coin: start in |0>, apply H, then the gates, then H. The chance of measuring 1.

    `gates_matrix` is the whole gate sequence multiplied into one matrix (see gate_product).
    """
    state = [1, 0]                             # start in |0>
    state = apply_gate(H, state)
    state = apply_gate(gates_matrix, state)
    state = apply_gate(H, state)
    amplitude_of_one = state[1]
    return abs(amplitude_of_one) ** 2          # the Born rule: chance = |amplitude| squared


def best_coins(target, max_t):
    """For each T budget, the coin whose chance of landing 1 is closest to `target`.

    Returns a list whose entry n is (chance of 1, gate list) for the best coin with at most
    n T gates.
    """

    def distance_from_target(matrix):
        # The score for a coin: how far its chance of 1 is from the target (0 is perfect).
        return abs(chance_of_one(matrix) - target)

    coins = []
    for gates in best_by_t_budget(distance_from_target, max_t):
        coins.append((chance_of_one(gate_product(gates)), gates))
    return coins


# ----- The search -----

# The 24 one-qubit Clifford gates, each written with H and S. Each one turns the Bloch
# sphere so that the x, y and z axes land on axes again, and there are 24 ways to do that
# (the 24 ways to set a cube down on a table). dial_demo_numpy.py finds this list by search.
CLIFFORD_SEQUENCES = [
    [],
    ['h'],
    ['s'],
    ['h', 's'],
    ['s', 'h'],
    ['s', 's'],
    ['h', 's', 'h'],
    ['h', 's', 's'],
    ['s', 'h', 's'],
    ['s', 's', 'h'],
    ['s', 's', 's'],
    ['h', 's', 'h', 's'],
    ['h', 's', 's', 'h'],
    ['h', 's', 's', 's'],
    ['s', 'h', 's', 's'],
    ['s', 's', 'h', 's'],
    ['h', 's', 'h', 's', 's'],
    ['h', 's', 's', 'h', 's'],
    ['s', 'h', 's', 's', 'h'],
    ['s', 'h', 's', 's', 's'],
    ['s', 's', 'h', 's', 's'],
    ['h', 's', 'h', 's', 's', 'h'],
    ['h', 's', 'h', 's', 's', 's'],
    ['h', 's', 's', 'h', 's', 's'],
]

# Each Clifford gate with its matrix, worked out once here because the search uses them
# over and over.
CLIFFORDS = []
for clifford_gates in CLIFFORD_SEQUENCES:
    CLIFFORDS.append((clifford_gates, gate_product(clifford_gates)))


def best_by_t_budget(score, max_t):
    """For each T budget 0, 1, ..., max_t: the gate list with the lowest score.

    `score` is a function that rates a gate sequence's matrix: lower is better. Returns a
    list whose entry n is the best gate list with at most n T gates.

    Trying every list of gates would take forever. A theorem saves us (Matsumoto and Amano,
    arXiv:0806.3834): every operator made of H, S and T can be written in exactly one
    standard form,

        a Clifford gate, then blocks of 'T H' or 'T H S', then maybe one more T,

    and this form has the fewest T gates possible. Each block holds one T gate. So we try
    each standard form once: every choice of blocks, with and without the last T, after
    each of the 24 Clifford gates.
    """
    best = {}                     # for each exact T count: (score, gate list), the best so far
    middles = [([], IDENTITY)]    # every choice of blocks so far, with its matrix; first, none
    for n_blocks in range(max_t + 1):
        for middle_gates, middle_matrix in middles:
            for last_gates in ([], ['t']):
                n_t = n_blocks + len(last_gates)
                if n_t > max_t:
                    continue
                # The last T acts after the blocks, so its matrix goes on the left.
                matrix = multiply(gate_product(last_gates), middle_matrix)
                found = best_after_a_clifford(score, middle_gates + last_gates, matrix)
                if n_t not in best or found[0] < best[n_t][0]:
                    best[n_t] = found
        if n_blocks < max_t:
            middles = one_block_longer(middles)
    return best_within_budget(best, max_t)


def best_after_a_clifford(score, gates, matrix):
    """Put each Clifford gate in front of a gate list, and return the best (score, gate list).

    `matrix` is the gate list's matrix. The Clifford gate acts first, so its matrix goes on
    the right of the product.
    """
    best_score = None
    best_gates = None
    for clifford_gates, clifford_matrix in CLIFFORDS:
        this_score = score(multiply(matrix, clifford_matrix))
        if best_score is None or this_score < best_score:
            best_score = this_score
            best_gates = clifford_gates + gates
    return best_score, best_gates


def one_block_longer(middles):
    """Every middle part made one block longer: once with 'T H' added, once with 'T H S'."""
    longer = []
    for middle_gates, middle_matrix in middles:
        for block in (['t', 'h'], ['t', 'h', 's']):
            # The block acts after the middle part, so its matrix goes on the left.
            longer.append((middle_gates + block, multiply(gate_product(block), middle_matrix)))
    return longer


def best_within_budget(best, max_t):
    """From the best gate list for each exact T count, the best within each budget.

    A budget of n T gates allows any count from 0 to n, so we keep the best seen so far.
    """
    within_budget = []
    record = None
    for n_t in range(max_t + 1):
        if record is None or best[n_t][0] < record[0]:
            record = best[n_t]
        within_budget.append(record[1])        # keep the gate list, drop its score
    return within_budget


# ----- The spin: turns about z -----

def turn_about_z(angle):
    """The matrix of a turn by `angle` radians about the z axis (Qiskit calls it Rz)."""
    return [[cmath.exp(-0.5j * angle), 0], [0, cmath.exp(0.5j * angle)]]


def rotation_error_degrees(gates_matrix, angle):
    """How far a gate sequence is from an exact turn of `angle` radians about z, in degrees.

    Undo the exact turn (turn back by the same angle) and see what is left over: for a
    perfect sequence, nothing. What is left is itself a turn, by some small angle about
    some axis, and that angle is the error. A turn by an angle 2a about the axis
    (nx, ny, nz) has the matrix

        | cos(a) - i nz sin(a)     (-i nx - ny) sin(a)  |
        | (-i nx + ny) sin(a)      cos(a) + i nz sin(a) |

    times an overall phase, which no measurement can see. Once the phase is divided out,
    sin(a) can be read off the top row. (Reading cos(a) off the diagonal fails for tiny
    errors: cos(a) is then so close to 1 that the computer rounds it to exactly 1.)
    """
    leftover = multiply(turn_about_z(-angle), gates_matrix)
    # Divide out the overall phase: the matrix above has determinant 1, so divide by the
    # square root of the leftover's determinant.
    determinant = leftover[0][0] * leftover[1][1] - leftover[0][1] * leftover[1][0]
    phase = cmath.sqrt(determinant)
    top_left = leftover[0][0] / phase          # cos(a) - i nz sin(a)
    top_right = leftover[0][1] / phase         # (-i nx - ny) sin(a)
    # |top right|^2 + (imaginary part of top left)^2 = (nx^2 + ny^2 + nz^2) sin(a)^2 = sin(a)^2
    sin_a = math.sqrt(abs(top_right) ** 2 + top_left.imag ** 2)
    sin_a = min(sin_a, 1.0)                    # rounding can push it a hair above 1
    return math.degrees(2 * math.asin(sin_a))


def best_rotations(angle, max_t):
    """For each T budget 0, 1, ..., max_t: the gate list closest to a turn of `angle` about z."""

    def error_of(matrix):
        # The score for a turn: its error in degrees (0 is perfect).
        return rotation_error_degrees(matrix, angle)

    return best_by_t_budget(error_of, max_t)


def qiskit_rotation(angle, epsilon):
    """The gate list Qiskit's compiler gives for a turn of `angle` radians about z.

    Qiskit uses the Ross-Selinger method, which finds a nearly shortest sequence with an
    error below `epsilon`. Qiskit 2.3.1 crashes on many angles when epsilon is larger than
    about 1e-4 (a bug in its gridsynth_rz). The lecture's 18-degree step is tested; if you
    try other angles, keep epsilon at 1e-6 or smaller.
    """
    circuit = gridsynth_rz(angle, epsilon)
    gates = []
    for instruction in circuit.data:
        gates.append(instruction.operation.name)
    return gates


def worst_drift(step, step_angle):
    """The spin's largest error over 20 time steps, computed exactly (no shots).

    After k steps, the chance of finding the spin pointing Left should be
    sin(k x step_angle / 2)^2. We follow the state step by step and compare at every step.
    """
    step_matrix = gate_product(step)
    state = apply_gate(H, [1, 0])              # the spin starts pointing Right: H|0>
    worst = 0
    for k in range(21):                        # after 0, 1, ..., 20 steps
        measured_state = apply_gate(H, state)  # H turns Right into |0> and Left into |1>
        chance_left = abs(measured_state[1]) ** 2
        exact_chance = math.sin(k * step_angle / 2) ** 2
        worst = max(worst, abs(chance_left - exact_chance))
        state = apply_gate(step_matrix, state)  # one more time step
    return worst


def plot_spin(step_angle, steps, chances_left):
    """Plot the spin's exact motion and, for each step in `steps`, its measured chances of Left.

    `chances_left[i]` holds the measured chances after 0, 1, ..., 20 steps of `steps[i]`.
    """
    fig, ax = plt.subplots(figsize=(9, 5.5), facecolor='#fcfcfb')
    ax.set_facecolor('#fcfcfb')
    smooth_k = []
    smooth_chance = []
    for i in range(400):
        k = 20 * i / 399
        smooth_k.append(k)
        smooth_chance.append(math.sin(k * step_angle / 2) ** 2)
    ax.plot(smooth_k, smooth_chance, color='#52514e', linewidth=1.5, zorder=3,
            label='the exact motion')
    colors = ['#2a78d6', '#eb6834', '#1baf7a']
    for i in range(len(steps)):
        ax.plot(range(21), chances_left[i], color=colors[i], linewidth=2, marker='o',
                markersize=8, markeredgecolor='#fcfcfb', markeredgewidth=2,
                label=f'{t_count(steps[i])} T gates per step')
    ax.set_title('A spin in a magnetic field, simulated with H, S and T gates', loc='left',
                 color='#0b0b0b')
    ax.set_xlabel('time step', color='#52514e')
    ax.set_ylabel('chance of finding the spin pointing Left', color='#52514e')
    ax.set_xticks(range(0, 21, 2))
    ax.set_ylim(-0.03, 1.03)
    ax.grid(color='#e1e0d9', linewidth=1)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color('#c3c2b7')
    ax.tick_params(colors='#898781')
    ax.legend(frameon=False, labelcolor='#0b0b0b')
    plt.show()
