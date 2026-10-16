
import numpy as np
import matplotlib.pyplot as plt
from integration import step_rk4
from integration import step_euler
from integration import step_milieu
from integration import step_heun
from integration import meth_n_step
from integration import meth_epsilon
from math import sin, cos
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp



def malthus_model(N, t):
    dN_dt = 2*np.array(N)           #Ici lambda = 2
    return dN_dt     

def verhulst_model(N, t):
    dN_dt = (2*np.array(N)*(1 - (np.array(N)/300))) #Ici lambda = 2, ĸ = 300, ce choix est cohérent pour des lapins
    return dN_dt
def Lotka_Volterra_model(y, t):
    # Ensure y is a NumPy array
    y = np.array(y)
    
    # Parameters for predator-prey model
    a, b = 1.2, 0.3  # Growth rate and interaction for prey
    c, d = 0.8, 1.0  # Growth rate and interaction for predators
    
    # Prey population (N): growth rate reduced by predation
    dN_dt = y[0] * (a - b * y[1])
    
    # Predator population (P): growth depends on prey availability
    dP_dt = y[1] * (c * y[0] - d)
    
    return np.array([dN_dt, dP_dt])
8
print("Malthus model with n_step and euler model : ")
print(meth_n_step([10], 0, 100, 0.1, malthus_model, step_euler))

print("Verhulst model with n_step and euler model : ")
print(meth_n_step([10], 0, 100, 0.1, verhulst_model, step_euler))



initial_conditions = [10, 5]
print("Lotka-Volterra model with n_step and Euler method:")
euler_result = meth_n_step(initial_conditions, 0, 100, 0.1, Lotka_Volterra_model, step_euler)
print("Final state (Euler):", euler_result)
    
print("\nLotka-Volterra model with n_step and RK4 method:")
rk4_result = meth_n_step(initial_conditions, 0, 100, 0.1, Lotka_Volterra_model, step_rk4)
print("Final state (RK4):", rk4_result)



# Simulation avec suivi de l'évolution
def meth_n_step_tracking(y0, t0, N, h, f, meth):
    y = np.array(y0, dtype=float)
    t = t0
    y_history = [y.copy()]
    t_history = [t]
    for _ in range(N):
        y = meth(y, t, h, f)
        t += h
        y_history.append(y.copy())
        t_history.append(t)
    return np.array(y_history), np.array(t_history)

# Modèle de Lotka-Volterra
def Lotka_Volterra_model(y, t):
    a, b = 1.2, 0.3  # Paramètres pour les proies
    c, d = 0.8, 1.0  # Paramètres pour les prédateurs
    dN_dt = y[0] * (a - b * y[1])       # Variation des proies
    dP_dt = y[1] * (c * y[0] - d)       # Variation des prédateurs
    return np.array([dN_dt, dP_dt])

# Paramètres de simulation
initial_conditions = [10, 5]  # 1000 proies, 5 prédateurs
t0 = 0
N = 1000
h = 0.1

# Exécution de la simulation
y_history, t_history = meth_n_step_tracking(initial_conditions, t0, N, h, Lotka_Volterra_model, step_rk4)

# Tracé des résultats
plt.figure(figsize=(10, 4))

# Évolution temporelle
plt.subplot(1, 2, 1)
plt.plot(t_history, y_history[:, 0], label='Proies (N)')
plt.plot(t_history, y_history[:, 1], label='Prédateurs (P)')
plt.title('Évolution des populations')
plt.xlabel('Temps')
plt.ylabel('Nombre d\'individus')
plt.legend()
plt.grid(True)

# Portrait de phase
plt.subplot(1, 2, 2)
plt.plot(y_history[:, 0], y_history[:, 1])
plt.title('Portrait de phase')
plt.xlabel('Proies (N)')
plt.ylabel('Prédateurs (P)')
plt.grid(True)

plt.tight_layout()
plt.show()
