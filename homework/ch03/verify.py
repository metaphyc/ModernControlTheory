"""Exact symbolic checks for the numerical feedback exercises in chapter 3."""
import sympy as sp

s = sp.symbols("s")


def charpoly(A):
    return A.charpoly(s).as_expr()


A = sp.Matrix([[0, 1, 0, 0], [3, 0, 0, 2], [0, 0, 0, 1], [0, -2, 0, 0]])
b = sp.Matrix([0, 0, 0, 1])
assert sp.Matrix.hstack(*(A**j * b for j in range(4))).det() == -12
K = sp.Matrix([[-16, -sp.Rational(22, 3), sp.Rational(5, 3), -6]])
assert sp.expand(charpoly(A + b*K) - (s+1)**2*((s+2)**2+1)) == 0
k1, k2, k3, k4 = sp.symbols("k1 k2 k3 k4")
generic_K = sp.Matrix([[k1, k2, k3, k4]])
cofactor_expansion = s*((s*s-3)*(s-k4)+2*s*(2-k2)-2*k1)-k3*(s*s-3)
assert sp.expand(charpoly(A+b*generic_K)-cofactor_expansion) == 0

A = sp.Matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
B = sp.Matrix([[0, 0], [0, 0], [1, 0], [0, 1]])
K = sp.Matrix([[-4, -6, -5, 0], [0, 0, 0, -3]])
assert sp.expand(charpoly(A + B*K) - (s+2)**2*((s+1)**2+1)) == 0

A = sp.Matrix([[0, 1, 0], [0, 0, 1], [8, 6, -3]])
b, c = sp.Matrix([0, 0, 1]), sp.Matrix([[6, 5, 1]])
K = sp.Matrix([[-24, -26, -5]])
go = (c*(s*sp.eye(3)-A).inv()*b)[0]
gc = (c*(s*sp.eye(3)-A-b*K).inv()*b)[0]
assert sp.cancel(go-(s+2)*(s+3)/((s+1)*(s-2)*(s+4))) == 0
assert sp.cancel(gc-(s+3)/((s+2)*(s+4))) == 0
Af = A+b*K
assert sp.Matrix.vstack(*(c*Af**j for j in range(3))).rank() == 2
unobservable_mode = sp.Matrix([1, -2, 4])
assert Af*unobservable_mode == -2*unobservable_mode
assert c*unobservable_mode == sp.zeros(1, 1)

A = sp.Matrix([[2, 1, 0], [0, 1, 0], [1, 0, 1]])
b, K = sp.Matrix([0, 1, 0]), sp.Matrix([[-36, -10, -24]])
T = sp.Matrix([[-4, -3, -2], [20, 12, 6], [1, 1, 1]])
assert T.det() == 2
assert (A+b*K)*T == T*sp.diag(-3, -2, -1)
generic_K = sp.Matrix([[k1, k2, k3]])
cofactor_expansion = (s-2)*(s-1-k2)*(s-1)-k1*(s-1)-k3
assert sp.expand(charpoly(A+b*generic_K)-cofactor_expansion) == 0
print("Chapter 3: all feedback, transfer-function, and similarity checks passed.")
