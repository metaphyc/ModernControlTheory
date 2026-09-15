"""Exact checks for the marked observer and output-feedback exercises.

Run: python3 homework/ch06/verify.py  (requires SymPy).
"""

import sympy as sp

s = sp.symbols("s")


def zero(expr):
    if isinstance(expr, sp.MatrixBase):
        assert all(sp.simplify(entry) == 0 for entry in expr), expr
    else:
        assert sp.simplify(expr) == 0, expr


def check_reduced_observer(A, B, L, expected_observer, K=None, expected_tf=None, implementation=None):
    """Check an observer against plant dynamics, and the actual controller interconnection."""
    n = A.rows
    C = sp.zeros(1, n)
    C[0, 0] = 1
    A11, A12 = A[:1, :1], A[:1, 1:]
    A21, A22 = A[1:, :1], A[1:, 1:]
    B1, B2 = B[:1, :], B[1:, :]
    F = A22 - L * A12
    H = F * L + A21 - L * A11
    E = B2 - L * B1
    # eta = z - L*y is the true counterpart of observer state w.
    T = (-L).row_join(sp.eye(n - 1))
    zero(T * A - F * T - H * C)
    zero(T * B - E)
    zero(F.charpoly(s).as_expr() - expected_observer)
    if K is None:
        return F, H, E
    Kz = K[:, 1:]
    Ky = K[:, :1] + Kz * L
    # Use the separately transcribed final controller coefficients, so a
    # correct observer formula cannot conceal an error in the displayed implementation.
    N, Hy, Ev, Cw, Dy = implementation
    zero(N - (F + E * Kz))
    zero(Hy - (H + E * Ky))
    zero(Ev - E)
    zero(Cw - Kz)
    zero(Dy - Ky)
    # Physical closed-loop state is [x;w], input is v.
    Ac = (A + B * Dy * C).row_join(B * Cw).col_join(
        (Hy * C).row_join(N)
    )
    Bc = B.col_join(Ev)
    Cc = C.row_join(sp.zeros(1, n - 1))
    # Transform to [x;e] with e = w - (z-L*y).
    S = sp.eye(n).row_join(sp.zeros(n, n - 1)).col_join(
        (-T).row_join(sp.eye(n - 1))
    )
    separated = (A + B * K).row_join(B * Kz).col_join(
        sp.zeros(n - 1, n).row_join(F)
    )
    zero(S * Ac * S.inv() - separated)
    zero(S * Bc - B.col_join(sp.zeros(n - 1, 1)))
    zero(Ac.charpoly(s).as_expr() - (A + B * K).charpoly(s).as_expr() * expected_observer)
    zero((Cc * (s * sp.eye(Ac.rows) - Ac).inv() * Bc)[0] - expected_tf)
    return F, H, E


# 1: full-order error dynamics.
A1 = sp.Matrix([[0, 1], [0, 0]])
C1 = sp.Matrix([[1, 0]])
L1 = sp.Matrix([6, 8])
zero((A1 - L1 * C1).charpoly(s).as_expr() - (s + 2) * (s + 4))
assert C1.col_join(C1 * A1).rank() == 2
print("6.1: observability and full-order observer poles verified")

# 2: reorder original state to [y=x2; z=x1].
A2 = sp.Matrix([[1, 2], [3, 1]])
B2 = sp.Matrix([2, 1])
F2, H2, E2 = check_reduced_observer(A2, B2, sp.Matrix([2]), s + 3)
zero(F2 - sp.Matrix([-3]))
zero(H2 - sp.Matrix([-5]))
zero(E2 - sp.Matrix([-3]))
print("6.2: implementable reduced observer and error dynamics verified")

# 6: full physical interconnection, including the two observer states.
A6 = sp.Matrix([[0, 1, 0], [0, 0, 1], [0, -2, -3]])
B6 = sp.Matrix([0, 0, 1])
C6 = sp.Matrix([[1, 0, 0]])
K6 = sp.Matrix([[-3, -2, -1]])
zero((C6 * (s * sp.eye(3) - A6).inv() * B6)[0] - 1 / (s * (s + 1) * (s + 2)))
assert sp.Matrix.hstack(B6, A6 * B6, A6**2 * B6).det() == -1
F6, H6, E6 = check_reduced_observer(
    A6, B6, sp.Matrix([7, 2]), (s + 5)**2, K6,
    1 / ((s + 3) * (s**2 + s + 1)),
    (sp.Matrix([[-7, 1], [-6, -4]]), sp.Matrix([-47, -53]),
     sp.Matrix([0, 1]), sp.Matrix([[-2, -1]]), sp.Matrix([[-19]])),
)
zero(F6 - sp.Matrix([[-7, 1], [-4, -3]]))
zero(H6 - sp.Matrix([-47, -34]))
zero((F6 + 5 * sp.eye(2))**2)
zero((A6 + B6 * K6).charpoly(s).as_expr() - (s + 3) * (s**2 + s + 1))
print("6.6: realization, repeated observer poles, separation and closed-loop transfer verified")

# 7: non-canonical plant and three-state observer.
A7 = sp.Matrix([[0, 1, 0, 0], [0, 0, -1, 0], [0, 0, 0, 1], [0, 0, 5, 0]])
B7 = sp.Matrix([0, 1, 0, -2])
C7 = sp.Matrix([[1, 0, 0, 0]])
K7 = sp.Matrix([[sp.Rational(4, 3), sp.Rational(10, 3), sp.Rational(49, 6), sp.Rational(25, 6)]])
assert sp.Matrix.hstack(*(A7**i * B7 for i in range(4))).det() == 36
assert sp.Matrix.vstack(*(C7 * A7**i for i in range(4))).det() == 1
pc = (s + 1) * (s + 2) * (s**2 + 2 * s + 2)
F7, H7, E7 = check_reduced_observer(
    A7, B7, sp.Matrix([9, -36, -84]), (s + 3) * (s**2 + 6 * s + 13), K7,
    (s**2 - 3) / pc,
    (sp.Matrix([[sp.Rational(-17, 3), sp.Rational(43, 6), sp.Rational(25, 6)],
                [36, 0, 1], [sp.Rational(232, 3), sp.Rational(-34, 3), sp.Rational(-25, 3)]]),
     sp.Matrix([sp.Rational(-1973, 3), 240, sp.Rational(5404, 3)]),
     sp.Matrix([1, 0, -2]),
     sp.Matrix([[sp.Rational(10, 3), sp.Rational(49, 6), sp.Rational(25, 6)]]),
     sp.Matrix([[sp.Rational(-1838, 3)]])),
)
zero(F7 - sp.Matrix([[-9, -1, 0], [36, 0, 1], [84, 5, 0]]))
zero(H7 - sp.Matrix([-45, 240, 576]))
zero((A7 + B7 * K7).charpoly(s).as_expr() - pc)
assert sp.gcd(s**2 - 3, pc) == 1
# Check the scalar input-output elimination displayed in the detailed solution.
k1, k2, k3, k4 = K7
eliminated = (s**2 - k2 * s - k1) * (s**2 - 3) - 2 * s**2 * (1 - k3 - k4 * s)
zero(eliminated - pc)
print("6.7: gains, separation, all seven poles and uncancelled transfer poles verified")
