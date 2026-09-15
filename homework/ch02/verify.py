"""Exact algebra checks for the computations in chapter 2.

Run with Python 3 and SymPy installed: python3 verify.py
The symbolic proofs in main.tex are not replaced by these checks.
"""
import sympy as sp


def controllability(a, b):
    return sp.Matrix.hstack(*(a**k * b for k in range(a.rows)))


def assert_zero(value):
    if isinstance(value, sp.MatrixBase):
        assert all(sp.simplify(entry) == 0 for entry in value), value
    else:
        assert sp.simplify(value) == 0, value


a, b, s = sp.symbols("a b s", real=True)
A = sp.Matrix([[-2, 1, 0], [0, -2, 0], [0, 0, -2]])
B = sp.Matrix([[a, 1], [2, 4], [b, 1]])
N = A + 2 * sp.eye(3)
assert N**2 == sp.zeros(3)
minor = sp.Matrix.hstack(B, (N * B)[:, 0])
assert_zero(minor.det() - 2 * (2 - 4 * b))
Q = controllability(A, B)
assert Q.subs(b, sp.Rational(1, 2)).rank() == 2
print("7(2): full-rank minor and exceptional controllability rank verified.")

A = sp.diag(0, 1)
B = sp.Matrix([1, 0])
T = sp.Matrix([[1, 1], [0, 1]])
Ah, Bh = T.inv() * A * T, T.inv() * B
assert Ah == sp.Matrix([[0, -1], [0, 1]])
assert_zero(controllability(Ah, Bh) - T.inv() * controllability(A, B))
v = sp.Matrix([0, 1])
assert_zero(controllability(A, B).T * v)
assert_zero(controllability(Ah, Bh).T * v)
assert controllability(Ah, Bh).T * (T.inv() * v) != sp.zeros(2, 1)
assert_zero(controllability(Ah, Bh).T * (T.T * v))
print("20: transformed reachable space and orthogonal-complement counterexample verified.")

A = sp.Matrix([[-1, -2, -2], [0, -1, 1], [1, 0, 1]])
B = sp.Matrix([2, 0, 1])
C = sp.Matrix([[0, 1, 1]])
T = sp.Matrix([[-6, -2, 2], [3, 1, 0], [3, 4, 1]])
Ac = sp.Matrix([[0, 1, 0], [0, 0, 1], [-3, -1, -1]])
Bc = sp.Matrix([0, 0, 1])
Cc = sp.Matrix([[6, 5, 1]])
assert controllability(A, B).det() == -18
assert T.det() == 18
assert_zero(A * T - T * Ac)
assert_zero(T * Bc - B)
assert_zero(C * T - Cc)
assert_zero((s * sp.eye(3) - A).det() - (s**3 + s**2 + s + 3))
original_tf = (C * (s * sp.eye(3) - A).inv() * B)[0]
canonical_tf = (Cc * (s * sp.eye(3) - Ac).inv() * Bc)[0]
assert_zero(original_tf - canonical_tf)
assert_zero(original_tf - (s**2 + 5 * s + 6) / (s**3 + s**2 + s + 3))
print("22: nonsingularity, coordinate transformation, characteristic polynomial and transfer function verified.")

# Check the endpoint-control construction of problem 15 on a nontrivial
# double integrator with a known, time-dependent forcing term.
t = sp.symbols("t", real=True)
A = sp.Matrix([[0, 1], [0, 0]])
B = sp.Matrix([0, 1])
transition = (A * (1 - t)).exp()
f = sp.Matrix([t, 1])
x0, xf = sp.Matrix([2, -1]), sp.Matrix([-1, 3])
W = sp.integrate(transition * B * B.T * transition.T, (t, 0, 1))
df = sp.integrate(transition * f, (t, 0, 1))
r = xf - A.exp() * x0 - df
control = (B.T * transition.T * W.inv() * r)[0]
endpoint = A.exp() * x0 + sp.integrate(transition * (B * control + f), (t, 0, 1))
assert W.det() > 0
assert_zero(endpoint - xf)
print("15: constructed control reaches a specified terminal state with time-dependent forcing.")
print("All chapter 2 algebra checks passed.")
