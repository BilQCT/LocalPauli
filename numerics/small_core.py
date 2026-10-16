# Independent enumeration implementation adopted for the numerical calculations.
"""Independent operator enumeration and exact-certificate helpers for n <= 3.

Coordinates are r_P = Tr(P rho), in lexicographic tensor order over I,X,Y,Z.
An integer column v therefore represents A = sum_P v_P P / 2**n.
Only numpy, scipy and sympy are required. No manuscript or legacy audit files
are read, and this module does not write files.
"""

from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations, product
from math import lcm, sqrt

import numpy as np
import sympy as sp
from scipy import sparse
from scipy.optimize import linprog


SQRT2 = sp.sqrt(2)


def pauli_data(n):
    """Binary XOR labels with 0=I, 1=X, 2=Y, 3=Z on each qubit."""
    digits = np.array(list(product(range(4), repeat=n)), dtype=np.int64)
    count = 4**n
    comm = np.zeros((count, count), dtype=np.uint8)
    phase = np.zeros((count, count), dtype=np.int8)
    # P_a P_b = i**phase[a,b] P_(a XOR b).
    local_phase = np.array([[0, 0, 0, 0], [0, 0, 1, 3],
                            [0, 3, 0, 1], [0, 1, 3, 0]])
    for j in range(n):
        a, b = digits[:, None, j], digits[None, :, j]
        comm ^= ((a != b) & (a != 0) & (b != 0)).astype(np.uint8)
        phase = (phase + local_phase[a, b]) % 4
    return digits, comm, phase


def isotropic_spaces(n, comm):
    """All isotropic binary subspaces, grouped by dimension, with one basis."""
    levels = [{(0,): ()}]
    for _ in range(n):
        next_level = {}
        for space, basis in levels[-1].items():
            candidates = (np.flatnonzero(~np.any(comm[list(basis)], axis=0))
                          if basis else range(4**n))
            for a in candidates:
                if a in space:
                    continue
                enlarged = tuple(sorted(space + tuple(int(a) ^ b for b in space)))
                next_level.setdefault(enlarged, basis + (int(a),))
        levels.append(next_level)
    return levels


def commuting_signs(basis, generator_signs, phase):
    signs = {0: 1}
    for a, sign in zip(basis, generator_signs):
        signs.update({a ^ b: int(s * sign * (1 if phase[a, b] == 0 else -1))
                      for b, s in list(signs.items())})
    return signs


def stabilizer_columns(n, levels, phase):
    rows, cols, values = [], [], []
    column = 0
    for _, basis in levels[n].items():
        for generator_signs in product([1, -1], repeat=n):
            signs = commuting_signs(basis, generator_signs, phase)
            rows.extend(signs)
            cols.extend([column] * len(signs))
            values.extend(signs.values())
            column += 1
    return sparse.csc_matrix((values, (rows, cols)),
                             shape=(4**n, column), dtype=np.int64)


def anticommuting_cliques(vertices, size, comm):
    def visit(prefix, remaining):
        if len(prefix) == size:
            yield prefix
            return
        if len(prefix) + len(remaining) < size:
            return
        for i, a in enumerate(remaining):
            yield from visit(prefix + (a,),
                             [b for b in remaining[i + 1:] if comm[a, b]])
    yield from visit((), vertices)


def cnc_columns(n, levels, comm, phase):
    """Enumerate every maximal CNC operator, not a sampled subset.

    A maximal CNC support is I union (a_1+I) union ... union (a_(2m+1)+I),
    where dim(I)=n-m and the nonzero quotient labels a_j in I^perp/I are
    pairwise anticommuting. All compatible value assignments are included.
    """
    rows, cols, values = [], [], []
    column = 0
    counts_by_m = {}
    for m in range(1, n + 1):
        starting_column = column
        for center, basis in levels[n - m].items():
            perp = (np.flatnonzero(~np.any(comm[list(basis)], axis=0))
                    if basis else range(4**n))
            representatives = sorted({min(int(a) ^ b for b in center)
                                      for a in perp} - {0})
            for clique in anticommuting_cliques(representatives, 2*m + 1, comm):
                for center_signs in product([1, -1], repeat=n - m):
                    signs = commuting_signs(basis, center_signs, phase)
                    for clique_signs in product([1, -1], repeat=2*m + 1):
                        rows.extend(signs)
                        cols.extend([column] * len(signs))
                        values.extend(signs.values())
                        for a, sign in zip(clique, clique_signs):
                            rows.extend(a ^ b for b in signs)
                            cols.extend([column] * len(signs))
                            values.extend(int(s * sign *
                                              (1 if phase[a, b] == 0 else -1))
                                          for b, s in signs.items())
                        column += 1
        counts_by_m[str(m)] = column - starting_column
    matrix = sparse.csc_matrix((values, (rows, cols)),
                                shape=(4**n, column), dtype=np.int64)
    return matrix, counts_by_m


def deterministic_columns(n):
    local = np.array([[1, *s] for s in product([1, -1], repeat=3)],
                     dtype=np.int64).T
    matrix = local
    for _ in range(n - 1):
        matrix = np.kron(matrix, local)
    return sparse.csc_matrix(matrix)


def single_qubit_clifford_actions():
    """Exactly the 24 proper signed permutations generated by H and S.

    These matrices act on the coefficient column (X,Y,Z). They do not
    include improper rotations or general outcome relabelings.
    """
    h = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], dtype=np.int8)
    s = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=np.int8)
    identity = np.eye(3, dtype=np.int8)
    found, queue = {identity.tobytes()}, [identity]
    while queue:
        current = queue.pop()
        for generator in (h, s):
            candidate = generator @ current
            key = candidate.tobytes()
            if key not in found:
                found.add(key)
                queue.append(candidate)
    actions = [np.frombuffer(key, dtype=np.int8).reshape(3, 3)
               for key in sorted(found)]
    assert len(actions) == 24
    assert all(round(np.linalg.det(a)) == 1 for a in actions)
    return actions


def gamma2_seed(n):
    """The paper's convention: each product is mod 3, the SUM is mod 2.

    X,Y,Z are identified with inputs 0,1,2. Reducing the entire sum mod 3
    before mod 2 would define a different function.
    """
    inputs = np.array(list(product(range(3), repeat=n)), dtype=np.int64)
    exponent = np.zeros(3**n, dtype=np.int64)
    for i, j in combinations(range(n), 2):
        exponent += (inputs[:, i] * inputs[:, j]) % 3
    return (1 - 2 * (exponent % 2)).astype(np.int8)


def maxweight_clifford_columns(n, digits):
    seed = gamma2_seed(n)
    actions = single_qubit_clifford_actions()
    orbit = set()
    for factors in product(actions, repeat=n):
        action = np.array([[1]], dtype=np.int8)
        for factor in factors:
            action = np.kron(action, factor)
        orbit.add((action @ seed).tobytes())
    full_weight = np.all(digits != 0, axis=1)
    matrix = np.zeros((4**n, len(orbit)), dtype=np.int64)
    matrix[0] = 1
    matrix[full_weight] = np.array([np.frombuffer(key, dtype=np.int8)
                                    for key in sorted(orbit)]).T
    return sparse.csc_matrix(matrix)


def resource_coordinates(n, edges, decorated, digits):
    bits = np.array(list(product(range(2), repeat=n)))
    ket = np.ones(2**n, dtype=complex) / sqrt(2**n)
    for i, j in edges:
        ket *= (-1.) ** (bits[:, i] * bits[:, j])
    h = np.array([[1, 1], [1, -1]]) / sqrt(2)
    for i in decorated:
        gate = np.eye(1)
        for j in range(n):
            gate = np.kron(gate, h if i == j else np.eye(2))
        ket = gate @ ket
    ket *= np.exp(1j * np.pi / 4 * np.sum(bits, axis=1))
    paulis = [np.eye(2), np.array([[0, 1], [1, 0]]),
              np.array([[0, -1j], [1j, 0]]), np.diag([1, -1])]
    coordinates = []
    for ds in digits:
        operator = np.eye(1)
        for d in ds:
            operator = np.kron(operator, paulis[d])
        coordinates.append(np.vdot(ket, operator @ ket).real)
    return np.array(coordinates)


def exact_resource_coordinates(n, edges, decorated, digits):
    """Build the actual state symbolically; no rationalization of float data."""
    bits = list(product(range(2), repeat=n))
    ket = sp.Matrix([(-1)**sum(row[i]*row[j] for i, j in edges) / sp.sqrt(2**n)
                     for row in bits])
    h = sp.Matrix([[1, 1], [1, -1]]) / SQRT2
    for i in decorated:
        ket = sp.kronecker_product(*(h if i == j else sp.eye(2)
                                     for j in range(n))) * ket
    t = (1 + sp.I) / SQRT2
    ket = sp.Matrix([sp.simplify(ket[k] * t**sum(row))
                     for k, row in enumerate(bits)])
    paulis = [sp.eye(2), sp.Matrix([[0, 1], [1, 0]]),
              sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1)]
    result = []
    for ds in digits:
        operator = sp.kronecker_product(*(paulis[int(d)] for d in ds))
        result.append(sp.simplify((ket.conjugate().T * operator * ket)[0]))
    return sp.Matrix(result)


@dataclass
class LPResult:
    objective: float
    primal: np.ndarray
    dual: np.ndarray
    primal_residual: float
    dual_violation: float
    primal_l1: float

    def summary(self):
        return {"value": self.objective,
                "primal_l1": self.primal_l1,
                "primal_residual": self.primal_residual,
                "dual_violation": self.dual_violation,
                "duality_gap": abs(self.primal_l1 - self.objective)}


def solve_robustness(matrix, target):
    """Dual max target.y subject |matrix.T y|<=1, with primal extraction."""
    constraints = sparse.vstack([matrix.T, -matrix.T]).tocsc()
    result = linprog(-target, A_ub=constraints,
                     b_ub=np.ones(2*matrix.shape[1]),
                     bounds=[(None, None)] * matrix.shape[0], method="highs")
    if not result.success:
        raise RuntimeError(result.message)
    multipliers = -result.ineqlin.marginals
    primal = multipliers[:matrix.shape[1]] - multipliers[matrix.shape[1]:]
    answer = LPResult(float(-result.fun), primal, result.x,
                      float(np.max(np.abs(matrix @ primal - target))),
                      float(max(0, np.max(np.abs(matrix.T @ result.x)) - 1)),
                      float(np.sum(np.abs(primal))))
    if max(answer.primal_residual, answer.dual_violation,
           abs(answer.primal_l1 - answer.objective)) > 1e-8:
        raise AssertionError(f"Inadequate numerical certificate: {answer.summary()}")
    return answer


def qpair(expression):
    expanded = sp.expand(sp.simplify(expression))
    b = sp.Rational(expanded.coeff(SQRT2))
    a = sp.Rational(sp.simplify(expanded - b*SQRT2))
    return Fraction(int(a.p), int(a.q)), Fraction(int(b.p), int(b.q))


def encode_qpair(expression):
    a, b = qpair(expression)
    return {"rational": [a.numerator, a.denominator],
            "sqrt2": [b.numerator, b.denominator]}


def decode_qpair(value):
    return Fraction(*value["rational"]), Fraction(*value["sqrt2"])


def pair_sign(pair):
    a, b = pair
    if not a:
        return (b > 0) - (b < 0)
    if not b or (a > 0) == (b > 0):
        return 1 if a > 0 else -1
    if a*a > 2*b*b:
        return 1 if a > 0 else -1
    return 1 if b > 0 else -1


def exact_primal_weights(integer_columns, denominator, exact_target):
    """Solve rational and sqrt(2) components together, avoiding expression swell."""
    rhs = sp.Matrix([[sp.Rational(x.numerator, x.denominator)
                      for x in qpair(v)] for v in exact_target])
    lhs = sp.Matrix(integer_columns.tolist()) / denominator
    solution, parameters = lhs.gauss_jordan_solve(rhs)
    if len(parameters):
        solution = solution.subs({symbol: 0 for symbol in parameters})
    assert lhs * solution == rhs
    return [sp.simplify(solution[i, 0] + solution[i, 1]*SQRT2)
            for i in range(solution.rows)]


def rational_dual_certificate(matrix, dual, exact_target):
    fractions = [Fraction(float(x)).limit_denominator(1000000) for x in dual]
    denominator = lcm(*(x.denominator for x in fractions))
    integers = [x.numerator * (denominator // x.denominator) for x in fractions]
    maximum = max_integer_dual_constraint(matrix, integers)
    if maximum > denominator:
        raise AssertionError("Rationalized dual fails an exact constraint")
    objective = sp.simplify(sum(sp.Rational(x.numerator, x.denominator)*v
                                for x, v in zip(fractions, exact_target)))
    return {"denominator": denominator, "numerators": integers,
            "maximum_absolute_constraint_numerator": maximum,
            "objective": encode_qpair(objective),
            "objective_expression": str(objective)}


def max_integer_dual_constraint(matrix, integers):
    if max(map(abs, integers), default=0) * matrix.shape[0] < 2**62:
        products = matrix.astype(np.int64).T @ np.array(integers, dtype=np.int64)
        return int(np.max(np.abs(products)))
    transposed = matrix.T.tocsr()
    maximum = 0
    for start, stop in zip(transposed.indptr[:-1], transposed.indptr[1:]):
        value = sum(int(s)*integers[int(i)] for s, i in
                    zip(transposed.data[start:stop], transposed.indices[start:stop]))
        maximum = max(maximum, abs(value))
    return maximum


def graph_state_twirl(n, edges, digits):
    """Integer sum of local-Clifford conjugations fixing T^n|G>.

    Generators g_i = ((X+Y)/sqrt(2))_i * product_(j in N(i)) Z_j.
    Include graph automorphisms as well. These actions preserve C and D.
    Return the sum matrix and the number of summands, not rounded floats.
    """
    r = np.array([[1, 0, 0, 0], [0, 0, 1, 0],
                  [0, 1, 0, 0], [0, 0, 0, -1]], dtype=np.int64)
    z, identity = np.diag([1, -1, -1, 1]), np.eye(4, dtype=np.int64)
    group = [np.eye(4**n, dtype=np.int64)]
    for i in range(n):
        neighbors = {b if a == i else a for a, b in edges if i in (a, b)}
        generator = np.eye(1, dtype=np.int64)
        for j in range(n):
            generator = np.kron(generator, r if j == i else
                                 z if j in neighbors else identity)
        group += [generator @ element for element in list(group)]
    edge_set = {tuple(sorted(edge)) for edge in edges}
    automorphisms = [p for p in permutations(range(n))
                     if {tuple(sorted((p[a], p[b]))) for a, b in edges} == edge_set]
    total = np.zeros((4**n, 4**n), dtype=np.int64)
    for p in automorphisms:
        action = np.zeros_like(total)
        for column, ds in enumerate(digits):
            row = sum(int(ds[p[j]]) * 4**(n - 1 - j) for j in range(n))
            action[row, column] = 1
        for element in group:
            total += action @ element
    return total, len(group)*len(automorphisms)


def validate_cnc_column(column, comm, phase):
    if column[0] != 1 or not set(column).issubset({-1, 0, 1}):
        return False
    support = np.flatnonzero(column)
    return all(column[a ^ b] == column[a]*column[b]*
               (1 if phase[a, b] == 0 else -1)
               for a in support for b in support if not comm[a, b])


def validate_deterministic_column(column, digits):
    n = digits.shape[1]
    local = [column[[4**(n - 1 - q)*a for a in (1, 2, 3)]] for q in range(n)]
    return column[0] == 1 and all(
        column[i] == np.prod([1 if a == 0 else local[q][a - 1]
                              for q, a in enumerate(ds)])
        for i, ds in enumerate(digits)) and all(abs(x) == 1 for x in column)
