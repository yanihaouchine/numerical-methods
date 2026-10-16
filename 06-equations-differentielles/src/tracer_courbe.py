import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
from utils import f_ref_2D, f_ref_2D_solution
from integration import *

t0 = 0
tf = 2 * np.pi
N = 100
h = (tf - t0) / N
y0 = np.array([1.0, 0.0])
ts = np.linspace(t0, tf, N+1)

methods = {'Euler': step_euler, 'Milieu': step_milieu, 'Heun': step_heun, 'RK4': step_rk4}

plt.figure(figsize=(8, 8))

theta = np.linspace(0, 2 * np.pi, 500)
plt.plot(np.cos(theta), np.sin(theta), 'k--', label='Solution exacte')

for name, method in methods.items():
    ys = [y0]
    y = y0.copy()
    t = t0
    for _ in range(N):
        y = method(y, t, h, f_ref_2D)
        ys.append(y)
        t += h
    ys = np.array(ys)
    plt.plot(ys[:, 0], ys[:, 1], label=name)

plt.gca().set_aspect('equal')
plt.title("Comparaison des méthodes d'intégration")
plt.xlabel("$y_0$")
plt.ylabel("$y_1$")
plt.legend()
plt.grid()
plt.savefig("./rapport/img/comparaison_methods.png")
plt.show()