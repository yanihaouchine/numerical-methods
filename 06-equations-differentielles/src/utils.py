import numpy as np
from typing import Union, Callable

class eq_diff:
    def __init__(self, t0: Union[int, float], y0: Union[int, float, np.ndarray], f: Callable[[Union[int, float, np.ndarray], Union[int, float]], Union[int, float, np.ndarray]]):
        self.t0 = t0
        self.y0 = y0
        self.f = f

def f_ref_1D(y, t):
    return y / (1 + t**2)

def f_ref_1D_solution(t):
    return np.arctan(t) + np.log(1 + t**2) / 2

def f_ref_2D(y: np.ndarray, t: np.ndarray):
    return np.array([-y[1], y[0]])

def f_ref_2D_solution(t):
    return np.array([np.cos(t), np.sin(t)])