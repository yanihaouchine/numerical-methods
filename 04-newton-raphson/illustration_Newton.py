import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from NewtonRaphson import Newton_Raphson



# Plot convergence speed with and without backtracking
def plot_convergence_speed():
    # Define the function f and its Jacobian J
    def f(U):
        x = U[0]
        return np.array([x**3 - 2*x + 2], float)

    def J(U):
        x = U[0]
        return np.array([[3*x**2 - 2]], float)

    # Compute the true solution
    true_solution = fsolve(lambda x: x**3 - 2*x + 2, x0=[0])[0]

    # Initial guess
    U0 = np.array([-2.0])
    # Number of iterations
    N = 50
    # Tolerance
    epsilon = 1e-6

    # Run Newton-Raphson without backtracking
    solution_no_bt, iterations_no_bt, _ = Newton_Raphson(f, J, U0, N, epsilon, backtracking=False)

    # Run Newton-Raphson with backtracking
    solution_bt, iterations_bt, _ = Newton_Raphson(f, J, U0, N, epsilon, backtracking=True)

    # Compute error at each iteration
    print(true_solution)
    print(solution_bt)
    print(solution_no_bt)
    errors_no_bt = np.abs(iterations_no_bt )
    errors_bt = np.abs(iterations_bt )

    # Plot the convergence speed
    plt.figure(figsize=(8, 6))
    plt.semilogy(range(len(errors_no_bt)), errors_no_bt, 'r-o', label='Without Backtracking')
    plt.semilogy(range(len(errors_bt)), errors_bt, 'b-o', label='With Backtracking')
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence Speed with and without Backtracking')
    plt.legend()
    plt.grid()
    plt.show()

plot_convergence_speed()
