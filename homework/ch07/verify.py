"""Verify exact CARE solutions, stability, costs and trajectories for marked LQR problems.

Run: python3 homework/ch07/verify.py  (requires SymPy).
"""

import sympy as sp

s = sp.symbols("s")
t = sp.symbols("t", real=True)


def zero(expr):
    if isinstance(expr, sp.MatrixBase):
        assert all(sp.simplify(entry) == 0 for entry in expr), expr
    else:
        assert sp.simplify(sp.radsimp(sp.expand(expr))) == 0, expr


def verify_lqr(A, B, Q, R, P, x0, expected_cost, expected_poly):
    zero(P - P.T)
    zero(A.T * P + P * A - P * B * B.T * P / R + Q)
    assert P[0, 0].is_positive
    assert sp.expand(P.det()).is_positive
    F = B.T * P / R
    Ac = A - B * F
    zero(Ac.charpoly(s).as_expr() - expected_poly)
    # Independent cost check: P is the closed-loop observability Gramian
    # of the stage cost Q+F.T*R*F, so its quadratic form integrates J.
    zero(Ac.T * P + P * Ac + Q + F.T * R * F)
    zero((x0.T * P * x0)[0] - expected_cost)
    x1, x2, u = sp.symbols("x1 x2 u", real=True)
    x = sp.Matrix([x1, x2])
    Vdot = (x.T * (A.T * P + P * A) * x)[0] + 2 * (x.T * P * B)[0] * u
    zero((x.T * Q * x)[0] + R * u**2 + Vdot - R * (u + (F * x)[0])**2)
    return F, Ac


A2 = sp.Matrix([[0, 1], [0, 0]])
B2 = sp.Matrix([0, 1])
Q2 = sp.Matrix([[2, 1], [1, 1]])
a2, b2 = sp.sqrt(2), sp.sqrt(1 + 2 * sp.sqrt(2))
P2 = sp.Matrix([[a2 * b2 - 1, a2], [a2, b2]])
J2 = (4 + sp.sqrt(2)) * b2 + 4 * sp.sqrt(2) - 1
F2, Ac2 = verify_lqr(A2, B2, Q2, sp.Integer(1), P2, sp.Matrix([1, 2]), J2, s**2 + b2 * s + a2)
zero(F2 - sp.Matrix([[a2, b2]]))
omega = sp.sqrt(2 * sp.sqrt(2) - 1) / 2
x1t = sp.exp(-b2 * t / 2) * (sp.cos(omega * t) + (2 + b2 / 2) / omega * sp.sin(omega * t))
zero(x1t.subs(t, 0) - 1)
zero(sp.diff(x1t, t).subs(t, 0) - 2)
zero(sp.diff(x1t, t, 2) + b2 * sp.diff(x1t, t) + a2 * x1t)
print(f"7.2: CARE, positive definiteness, stable poles, time trajectory and J*={sp.N(J2, 10)} verified")

A5 = sp.diag(1, 2)
B5 = sp.Matrix([1, 1])
C5 = sp.Matrix([[1, 2]])
Q5 = C5.T * C5
a5, b5 = 2 * sp.sqrt(2), sp.sqrt(6) / 2
sigma, rho = a5 + b5, a5 * b5
f1, f2 = -1 - sigma - rho, 4 + 2 * sigma + rho
p = 22 + 10 * sp.sqrt(2) + 8 * sp.sqrt(3) + 9 * sp.sqrt(6)
q = -24 - 14 * sp.sqrt(2) - 12 * sp.sqrt(3) - 10 * sp.sqrt(6)
r = 32 + 22 * sp.sqrt(2) + 16 * sp.sqrt(3) + 12 * sp.sqrt(6)
P5 = sp.Matrix([[p, q], [q, r]])
J5 = 24 + 6 * sp.sqrt(2) + 8 * sp.sqrt(6)
x05 = sp.Matrix([2, 1])
F5, Ac5 = verify_lqr(A5, B5, Q5, sp.Integer(2), P5, x05, J5, (s + a5) * (s + b5))
assert B5.row_join(A5 * B5).det() == 1
assert C5.col_join(C5 * A5).det() == 2
zero(F5 - sp.Matrix([[f1, f2]]))
zero(p - (f1**2 - sp.Rational(1, 2)))
zero(q - (2 * f1 * f2 - 2) / 3)
zero(r - (f2**2 / 2 - 1))
H = A5.row_join(-B5 * B5.T / 2).col_join((-Q5).row_join(-A5.T))
zero(H.charpoly(s).as_expr() - (s**2 - 8) * (s**2 - sp.Rational(3, 2)))
# Check the displayed Schur-complement / rank-one determinant derivation.
D, E = s * sp.eye(2) - A5, s * sp.eye(2) + A5.T
rank_one_det = D.det() * E.det() * (
    1 - (C5 * D.inv() * B5)[0] * (B5.T * E.inv() * C5.T)[0] / 2
)
zero(rank_one_det - H.charpoly(s).as_expr())
zero(H.charpoly(s).as_expr() - (s**2 - 1) * (s**2 - 4) + (3*s - 4)*(3*s + 4)/2)
zero(H * sp.eye(2).col_join(P5) - sp.eye(2).col_join(P5) * Ac5)
ca = (rho * (3 + a5) + 2 * b5) / (a5 - b5)
cb = -(rho * (3 + b5) + 2 * a5) / (a5 - b5)
ut = ca * sp.exp(-a5 * t) + cb * sp.exp(-b5 * t)
zero(ut.subs(t, 0) + (F5 * x05)[0])
zero(sp.diff(ut, t).subs(t, 0) + (F5 * Ac5 * x05)[0])
zero(sp.diff(ut, t, 2) + sigma * sp.diff(ut, t) + rho * ut)
print(f"7.5: CARE, Hamiltonian stable subspace, exact feedback, time control and J*={sp.N(J5, 10)} verified")
