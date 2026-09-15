"""独立核对习题一中的传递函数、状态响应和离散化。依赖 SymPy。"""
import sympy as sp

s, z, t, tau = sp.symbols("s z t tau", real=True)


def zero(expr):
    if isinstance(expr, sp.MatrixBase):
        assert all(sp.simplify(v) == 0 for v in expr), expr
    else:
        assert sp.simplify(expr) == 0, expr


# 题 2：从模型恢复串联 RLC 的传递函数。
R, L, C = sp.symbols("R L C", positive=True)
A = sp.Matrix([[-R/L, -1/L], [1/C, 0]])
B = sp.Matrix([1/L, 0])
zero((sp.Matrix([[0, 1]]) * (s*sp.eye(2)-A).inv()*B)[0]
     - 1/(L*C*s**2+R*C*s+1))

# 题 4：零状态传递函数；保留的二阶模型还能编码全部初值。
A = sp.Matrix([[0, 1], [1, 0]])
B = sp.Matrix([1, 1])
zero((sp.Matrix([[1, 0]])*(s*sp.eye(2)-A).inv()*B)[0] - 1/(s-1))
assert sp.Matrix.hstack(B, A*B).rank() == 1
assert sp.Matrix.vstack(sp.Matrix([[1, 0]]), sp.Matrix([[1, 0]])*A).rank() == 2

# 题 7：使用原始矩阵直接计算，核对正文给出的各个传函分量。
A = sp.Matrix([[0, 1, 0], [-2, -3, 0], [-1, 1, 3]])
B = sp.Matrix([0, 1, 2])
zero((sp.Matrix([[0, 0, 1]])*(s*sp.eye(3)-A).inv()*B)[0]
     - (2*s**2+7*s+3)/((s-3)*(s+1)*(s+2)))
A = sp.Matrix([[0, 1, 0], [0, 0, 1], [-3, -1, -2]])
B = sp.Matrix([[1, 0], [0, 1], [1, 1]])
zero(sp.Matrix([[1, 1, 1]])*(s*sp.eye(3)-A).inv()*B
     - sp.Matrix([[2*s**2-1, 2*s**2+3*s]])/(s**3+2*s**2+s+3))

# 题 10：同时核对两个微分方程和初始条件。
A = sp.Matrix([[0, 1], [-2, -3]])
B = sp.Matrix([2, 0])
x = sp.Matrix([(4*t-1)*sp.exp(-t)+sp.exp(-2*t),
               (3-4*t)*sp.exp(-t)-2*sp.exp(-2*t)])
zero(x.diff(t)-A*x-B*sp.exp(-t))
zero(x.subs(t, 0)-sp.Matrix([0, 1]))

# 题 15：恢复原脉冲传递函数及一步递推的输入输出恒等式。
B = sp.Matrix([2, -3])
Cy = sp.Matrix([[1, 0]])
zero((Cy*(z*sp.eye(2)-A).inv()*B)[0] - (2*z+3)/(z**2+3*z+2))
x0 = sp.Matrix(sp.symbols("x1 x2"))
u0, u1 = sp.symbols("u0 u1")
x1 = A*x0+B*u0
x2 = A*x1+B*u1
zero((Cy*x2+3*Cy*x1+2*Cy*x0)[0]-2*u1-3*u0)

# 题 16：直接核对矩阵指数与保持器积分。
A = sp.Matrix([[0, 1], [0, 0]])
B = sp.Matrix([0, 1])
zero((2*A).exp()-sp.Matrix([[1, 2], [0, 1]]))
zero(sp.integrate((tau*A).exp()*B, (tau, 0, 2))-sp.Matrix([2, 2]))
print("ch01: symbolic checks passed")
