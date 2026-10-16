import numpy as np
import matplotlib.pyplot as plt
from integration import step_rk4
from math import sin, cos
g = 9.81  # gravitational constant

def simple_pendulum(L):
    def f(y, _):
        theta, omega = y
        return np.array([omega, -(g / L) * sin(theta)])
    return f

def double_pendulum(m1=1.0, m2=1.0, l1=1.0, l2=1.0):
    def f(y, _):
        theta1, omega1, theta2, omega2 = y
        delta = theta2 - theta1

        denom = (m1 + m2 * sin(delta)**2)

        domega1 = (
            m2 * l1 * omega1**2 * sin(delta) * cos(delta)
            + m2 * g * sin(theta2) * cos(delta)
            + m2 * l2 * omega2**2 * sin(delta)
            - (m1 + m2) * g * sin(theta1)
        ) / (l1 * denom)

        domega2 = (
            -m2 * l2 * omega2**2 * sin(delta) * cos(delta)
            + (m1 + m2) * (g * sin(theta1) * cos(delta) - l1 * omega1**2 * sin(delta) - g * sin(theta2))
        ) / (l2 * denom)

        return np.array([omega1, domega1, omega2, domega2])
    
    return f

def compute_period(theta0, L, t_max=20, method=step_rk4, h=0.01):
    y0 = np.array([theta0, 0])
    t0 = 0
    N = int(t_max / h)
    y = y0
    t = t0
    ys = [y]
    ts = [t]
    f = simple_pendulum(L)
    for _ in range(1, N):
        y = method(y, t, h, f)
        t += h
        ys.append(y)
        ts.append(t)

    thetas = [y[0] for y in ys]
    peaks = []
    for i in range(1, len(thetas)-1):
        if thetas[i-1] < thetas[i] > thetas[i+1]:
            peaks.append(ts[i])

    if len(peaks) >= 2:
        period = np.mean(np.diff(peaks))
        return period
    else:
        return None
    
def plot_frequency(L, theta_max, num_points=100):
    thetas = np.linspace(0.01, theta_max, num_points)
    periods = np.array([compute_period(theta, L) for theta in thetas])
    frequencies = 1 / periods
    
    plt.plot(thetas, frequencies)
    plt.xlabel('Angle initial θ0 (radians)')
    plt.ylabel('Fréquence d\'oscillation (Hz)')
    #plt.title(f'Fréquence d\'oscillation en fonction de θ0')
    freq_small_oscillations = 1/(2 * np.pi) * np.sqrt(g/L)
    plt.ylim(bottom=0)
    plt.axhline(y=freq_small_oscillations, color='r', linestyle='--', 
                label=f'Fréquence théorique sqrt(g/L) petites oscillations: {freq_small_oscillations:.4f} Hz')
    plt.grid(True)
    plt.legend()
    plt.savefig('rapport/img/pendulum_frequency')
    plt.show()

L = 1.0
theta_max = np.pi / 2

plot_frequency(L, theta_max)


def plot_extremity_double_pendulum():
    method=step_rk4
    l1 = l2 = 1.0
    y0 = np.array([np.pi / 2, 0.0, np.pi / 2 + 0.01, 0.0])
    t0 = 0.0
    h = 0.01
    N = 5000

    f = double_pendulum(l1=l1, l2=l2)

    ys = [y0]
    ts = [t0]
    y = y0
    t = t0

    for _ in range(1, N):
        y = method(y, t, h, f)
        t += h
        ys.append(y)
        ts.append(t)

    ys = np.array(ys)

    theta1 = ys[:, 0]
    theta2 = ys[:, 2]

    x1 = l1 * np.sin(theta1)
    y1 = -l1 * np.cos(theta1)

    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)

    plt.figure(figsize=(8, 6))
    plt.plot(x2, y2, color='blue', linewidth=1)
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    #plt.title("Trajectoire de l'extrémité du pendule à deux maillons")
    plt.grid(True)
    plt.axis("equal")
    plt.savefig('rapport/img/extremity_double_pendulum')
    plt.show()

plot_extremity_double_pendulum()

def compute_flip_time(theta1_0, theta2_0, l1=1.0, l2=1.0, t_max=20, h=0.4):
    y = np.array([theta1_0, 0.0, theta2_0, 0.0])
    f = double_pendulum(l1=l1, l2=l2)
    t = 0.0
    while t < t_max:
        y = step_rk4(y, t, h, f)
        t += h
        theta1 = y[0]
        if abs(theta1) > np.pi: # retournement
            return t
    return t_max

def flip_time_map(theta_range=np.pi, N=200, t_max=20):
    theta_vals = np.linspace(-theta_range, theta_range, N)
    flip_times = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            t_flip = compute_flip_time(theta_vals[i], theta_vals[j], t_max=t_max)
            flip_times[j, i] = t_flip

    plt.figure(figsize=(8, 6))
    plt.imshow(flip_times, extent=[-theta_range, theta_range, -theta_range, theta_range],
               origin='lower', cmap='inferno', aspect='auto')
    plt.colorbar(label="Temps de retournement (s)")
    plt.xlabel("θ₁ initial (rad)")
    plt.ylabel("θ₂ initial (rad)")
    #plt.title("Temps de premier retournement en fonction des angles initiaux")
    plt.savefig("rapport/img/first_flip_time")
    plt.show()

flip_time_map()
