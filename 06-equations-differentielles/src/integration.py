import numpy as np
from utils import eq_diff

def step_euler(y, t, h, f):
    return y + h * f(y, t)

def step_milieu(y, t, h, f):
    p1 = f(y + (h / 2) * f(y, t), t + (h / 2))
    return y + h * p1

def step_heun(y, t, h, f):
    p1 = f(y, t)
    p2 = f(y + (h * p1), t + h)
    return y + (h / 2) * (p1 + p2)

def step_rk4(y, t, h, f):
    p1 = f(y, t)
    p2 = f(y + (h / 2) * p1, t + (h / 2))
    p3 = f(y + (h / 2) * p2, t + (h / 2))
    p4 = f(y + h * p3, t + h)
    return y + (h / 6) * (p1 + 2 * p2 + 2 * p3 + p4)

def meth_n_step(y0, t0, N, h, f, meth):
    y = y0.copy()
    t = t0
    for i in range(N):
        y = meth(y, t, h, f)
        t += h
    return y

def meth_epsilon(y0, t0, tf, eps, f, meth):
    h = (tf - t0) / 100
    y = y0.copy()
    t = t0

    while t < tf:
        if t + h > tf:
            h = tf - t
        y1 = meth(y, t, h, f)
        y2_half = meth(meth(y, t, h/2, f), t + h/2, h/2, f)
        err = np.linalg.norm(y2_half - y1)

        if err < eps:
            t += h
            y = y2_half
            h *= 1.5
        else:
            h /= 2
    return t, y
