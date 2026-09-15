"""Check all selected chapter 5 realizations using exact rational arithmetic."""
import sympy as sp

s = sp.symbols("s")
R = sp.Rational


def transfer(A, B, C, D=None):
    if D is None:
        D = sp.zeros(C.rows, B.cols)
    return (C*(s*sp.eye(A.rows)-A).inv()*B+D).applyfunc(sp.cancel)


def ranks(A, B, C):
    n = A.rows
    controllability = sp.Matrix.hstack(*(A**j*B for j in range(n)))
    observability = sp.Matrix.vstack(*(C*A**j for j in range(n)))
    return controllability.rank(), observability.rank()


def check(A, B, C, target, D=None, expected=None):
    expected = expected or (A.rows, A.rows)
    assert ranks(A, B, C) == expected
    assert (transfer(A, B, C, D)-sp.Matrix(target)).applyfunc(sp.cancel).is_zero_matrix


# 1(c): full-state analysis and the actual coordinate change x = T z.
A = sp.Matrix([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, -1, 0], [0, 0, 0, 1]])
B = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
C = sp.Matrix([[1, 0, 0, 2], [0, 0, 2, 1]])
G = sp.Matrix([[1/(s*(s+1)), 2/(s-1)], [2/(s+1), 1/(s-1)]])
check(A, B, C, G, expected=(4, 3))
assert sp.Matrix.hstack(B, A*B).det() == 1
T = sp.eye(4)[:, [0, 2, 3, 1]]
At = sp.Matrix([[0, 1, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0]])
Bt = sp.Matrix([[0, 0], [1, 0], [0, 1], [0, 0]])
Ct = sp.Matrix([[1, 0, 2, 0], [0, 2, 1, 0]])
assert A*T == T*At and B == T*Bt and C*T == Ct
check(At[:3, :3], Bt[:3, :], Ct[:, :3], G)

# 4: direct term and the strict-proper remainders both matter.
A = sp.Matrix([[0, 1, 0], [0, 0, 1], [-48, -44, -12]])
B, C = sp.Matrix([0, 0, 1]), sp.Matrix([[3, 4, 1]])
check(A, B, C, [[(s+1)*(s+3)/((s+2)*(s+4)*(s+6))]])
A = sp.Matrix([[0, 1, 0], [0, 0, 1], [-R(3, 2), -3, -2]])
C = sp.Matrix([[-R(13, 4), -R(11, 2), -R(9, 2)]])
check(A, B, C, [[(5*s**3+s*s+4*s+1)/(2*s**3+4*s*s+6*s+3)]], sp.Matrix([[R(5, 2)]]))
observability = sp.Matrix.vstack(C, C*A, C*A*A)
assert observability == sp.Matrix([[-13, -22, -18], [27, 41, 14], [-21, -15, 13]])/4
assert observability.det() == -R(3677, 64)

# 5: independently check each controllable, observable, and Jordan realization.
Ga = (3*s*s+17*s+25)/((s+2)**3*(s+3))
Ac = sp.Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [-24, -44, -30, -9]])
Bc, Cc = sp.Matrix([0, 0, 0, 1]), sp.Matrix([[25, 17, 3, 0]])
Ao = sp.Matrix([[0, 0, 0, -24], [1, 0, 0, -44], [0, 1, 0, -30], [0, 0, 1, -9]])
Bo, Co = sp.Matrix([25, 17, 3, 0]), sp.Matrix([[0, 0, 0, 1]])
AJ = sp.Matrix([[-2, 1, 0, 0], [0, -2, 1, 0], [0, 0, -2, 0], [0, 0, 0, -3]])
BJ, CJ = sp.Matrix([0, 0, 1, 1]), sp.Matrix([[3, 2, 1, -1]])
for A, B, C in ((Ac, Bc, Cc), (Ao, Bo, Co), (AJ, BJ, CJ)):
    check(A, B, C, [[Ga]])
TJ = sp.Matrix([[1, -1, 1, -1], [-2, 3, -3, 3], [4, -8, 9, -9], [-8, 20, -26, 27]])
assert TJ.det() == 1
assert Ac*TJ == TJ*AJ and TJ*BJ == Bc and Cc*TJ == CJ
r = sp.symbols("r")
assert sp.expand(3*r*r+5*r+3-(r*r*(r+1)+2*r*(r+1)+3*(r+1)-r**3)) == 0

Gb = (s*s+8*s+15)/(s**3+7*s*s+14*s+8)
Ac = sp.Matrix([[0, 1, 0], [0, 0, 1], [-8, -14, -7]])
Bc, Cc = sp.Matrix([0, 0, 1]), sp.Matrix([[15, 8, 1]])
Ao = sp.Matrix([[0, 0, -8], [1, 0, -14], [0, 1, -7]])
Bo, Co = sp.Matrix([15, 8, 1]), sp.Matrix([[0, 0, 1]])
AJ, BJ = sp.diag(-1, -2, -4), sp.ones(3, 1)
CJ = sp.Matrix([[R(8, 3), -R(3, 2), -R(1, 6)]])
for A, B, C in ((Ac, Bc, Cc), (Ao, Bo, Co), (AJ, BJ, CJ)):
    check(A, B, C, [[Gb]])
TJ = sp.Matrix([[R(1, 3), -R(1, 2), R(1, 6)], [-R(1, 3), 1, -R(2, 3)], [R(1, 3), -2, R(8, 3)]])
assert TJ.det() == R(1, 6)
assert Ac*TJ == TJ*AJ and TJ*BJ == Bc and Cc*TJ == CJ

# 6(4): distinct scalar poles require three states.
A = sp.Matrix([[-1, 0, 0], [0, 0, 1], [0, -2, -2]])
B, C = sp.Matrix([[1, 0], [0, 0], [0, 1]]), sp.Matrix([[1, 1, 0]])
check(A, B, C, [[1/(s+1), 1/(s*s+2*s+2)]])
assert sp.Matrix.vstack(C, C*A, C*A*A) == sp.Matrix([[1, 1, 0], [-1, 0, 1], [1, -2, -2]])

# 6(7): residue ranks prove the minimal dimension before checking a realization.
G = sp.Matrix([[1/(s*(s+1)), 2/(s+2)], [2/(s+1), 1/(s+1)]])
residues = [G.applyfunc(lambda g: sp.residue(g, s, pole)) for pole in (0, -1, -2)]
assert [residue.rank() for residue in residues] == [1, 2, 1]
A = sp.diag(0, -1, -1, -2)
B = sp.Matrix([[1, 0], [1, 0], [0, 1], [0, 1]])
C = sp.Matrix([[1, -1, 0, 2], [0, 2, 1, 0]])
check(A, B, C, G)

# 10: realization, positive feedback convention, and closed-loop transfer function.
A, B, C = sp.Matrix([[0, 1], [0, -5]]), sp.Matrix([0, 1]), sp.Matrix([[100, 0]])
check(A, B, C, [[100/(s*(s+5))]])
K = sp.Matrix([[-800, -35]])
Af = A+B*K
assert sp.expand(Af.charpoly(s).as_expr()-(s*s+40*s+800)) == 0
check(Af, B, C, [[100/(s*s+40*s+800)]])
print("Chapter 5: all realizations, ranks, residues, transformations, and poles passed.")
