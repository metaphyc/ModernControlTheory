"""Exact algebra checks for chapter 4. Requires Python 3 and SymPy."""
import sympy as sp


def controllability(a, b):
    return sp.Matrix.hstack(*(a**k * b for k in range(a.rows)))


def observability(a, c):
    return sp.Matrix.vstack(*(c * a**k for k in range(a.rows)))


def assert_zero(value):
    if isinstance(value, sp.MatrixBase):
        assert all(sp.simplify(entry) == 0 for entry in value), value
    else:
        assert sp.simplify(value) == 0, value


a, b, s, k = sp.symbols("a b s k", real=True)
A = sp.Matrix([[-1, 1, a], [0, -2, 1], [0, 0, -3]])
B = sp.Matrix([0, 0, 1])
C = sp.Matrix([[0, 0, 1]])
assert_zero(controllability(A, B).det() + a + 1)
assert observability(A, C).rank() == 1
A = sp.Matrix([[0, 0, 1], [0, 1, 0], [-2, -3, -5]])
B = sp.Matrix([0, 1, a])
C = sp.Matrix([[0, 1, b]])
Q, O = controllability(A, B), observability(A, C)
assert_zero(Q.det() - (8 * a**2 + 21 * a + 9))
assert_zero(O.det() - 2 * b**2 * (3 * b - 8))
for special_a in [(-21 + 3 * sp.sqrt(17)) / 16, (-21 - 3 * sp.sqrt(17)) / 16]:
    assert_zero(Q.det().subs(a, special_a))
assert O.subs(b, 0).rank() == 1
assert O.subs(b, sp.Rational(8, 3)).rank() == 2
assert Q.subs(a, 0).rank() == O.subs(b, 1).rank() == 3
print("6: symbolic determinant conditions and exceptional observability ranks verified.")

A = sp.Matrix([[-1, -2, -2], [0, -1, 1], [1, 0, 1]])
B = sp.Matrix([2, 0, 1])
C = sp.Matrix([[1, 1, 0]])
S = sp.Matrix([[0, 3, -3], [0, -2, -1], [1, 1, 0]])
T = sp.Matrix([[-sp.Rational(1, 9), sp.Rational(1, 3), 1],
               [sp.Rational(1, 9), -sp.Rational(1, 3), 0],
               [-sp.Rational(2, 9), -sp.Rational(1, 3), 0]])
Ao = sp.Matrix([[0, 0, -3], [1, 0, -1], [0, 1, -1]])
Bo = sp.Matrix([-3, -1, 2])
Co = sp.Matrix([[0, 0, 1]])
assert observability(A, C).det() == 9
assert S.det() == -9
assert_zero(S * T - sp.eye(3))
assert_zero(S * A - Ao * S)
assert_zero(S * B - Bo)
assert_zero(C * T - Co)
original_tf = (C * (s * sp.eye(3) - A).inv() * B)[0]
canonical_tf = (Co * (s * sp.eye(3) - Ao).inv() * Bo)[0]
assert_zero(original_tf - canonical_tf)
assert_zero(original_tf - (2 * s**2 - s - 3) / (s**3 + s**2 + s + 3))
print("13: observability, inverse transformation and transfer-function equivalence verified.")

# Check both PBH factorization identities from problem 17 with symbolic
# rectangular B, C and F, rather than merely sampling feedback gains.
A = sp.Matrix(2, 2, sp.symbols("a0:4"))
B = sp.Matrix(2, 2, sp.symbols("b0:4"))
C = sp.Matrix(1, 2, sp.symbols("c0:2"))
F = sp.Matrix(2, 1, sp.symbols("f0:2"))
Acl = A + B * F * C
right_multiplier = sp.BlockMatrix([[sp.eye(2), sp.zeros(2)], [-F * C, sp.eye(2)]]).as_explicit()
left_multiplier = sp.BlockMatrix([[sp.eye(2), -B * F], [sp.zeros(1, 2), sp.eye(1)]]).as_explicit()
assert_zero(sp.Matrix.hstack(s * sp.eye(2) - A, B) * right_multiplier
            - sp.Matrix.hstack(s * sp.eye(2) - Acl, B))
assert_zero(left_multiplier * sp.Matrix.vstack(s * sp.eye(2) - A, C)
            - sp.Matrix.vstack(s * sp.eye(2) - Acl, C))
assert right_multiplier.det() == left_multiplier.det() == 1
print("17: both symbolic PBH factorizations and multiplier invertibility verified.")

A = sp.Matrix([[0, 1], [-1, 0]])
B = sp.Matrix([0, 1])
C = sp.Matrix([[1, 0]])
Acl = A + B * k * C
assert controllability(A, B).rank() == observability(A, C).rank() == 2
assert Acl.trace() == 0
assert_zero((s * sp.eye(2) - Acl).det() - (s**2 + 1 - k))
assert Acl.subs(k, 1) != sp.zeros(2)
assert Acl.subs(k, 1)**2 == sp.zeros(2)
full_state_feedback = A + B * sp.Matrix([[-1, -2]])
assert_zero((s * sp.eye(2) - full_state_feedback).det() - (s**2 + 2 * s + 2))
assert set(full_state_feedback.eigenvals()) == {-1 + sp.I, -1 - sp.I}
print("18: open-loop ranks, closed-loop characteristic polynomial and k=1 Jordan degeneracy verified.")
print("18: the illustrative full-state feedback has poles -1 +/- i.")
print("All chapter 4 algebra checks passed.")
