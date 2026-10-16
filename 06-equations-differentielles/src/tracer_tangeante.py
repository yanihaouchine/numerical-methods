import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib import pyplot as plt
from utils import *
import numpy as np

N = 20

x_range=(-5, 5)
y_range=(-5, 5)
f = f_ref_2D

x = np.linspace(x_range[0], x_range[1], N)
y = np.linspace(y_range[0], y_range[1], N)
X, Y = np.meshgrid(x, y)

U, V = np.zeros_like(X), np.zeros_like(Y)
for i in range(N):
    for j in range(N):
        y_vec = np.array([X[i, j], Y[i, j]])
        t_vec = np.array([X[i, j], Y[i, j]])
        dy_dt = f(y_vec, t_vec)
        U[i, j] = dy_dt[0]
        V[i, j] = dy_dt[1]

magnitude = np.sqrt(U**2 + V**2)
U /= magnitude
V /= magnitude

plt.gca().set_aspect('equal')
plt.quiver(X, Y, U, V, angles='xy')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Champ des tangentes')
plt.xlim(x_range)
plt.ylim(y_range)
plt.savefig('./rapport/img/champ_tangente.png')
plt.show()
