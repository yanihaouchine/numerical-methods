import numpy as np

def Newton_Raphson(f, J, U0, N, epsilon, backtracking=True):
    U = U0
    f_u = [np.linalg.norm(f(U))]
    for i in range(N):
        if np.linalg.norm(f(U)) < epsilon:
            break
        delta, _, _, _ = np.linalg.lstsq(J(U), -f(U), rcond=None)
        U_new = U + delta
        if backtracking:
            # Backtracking
            alpha = 1.0
            rho = 0.5
            norme = np.linalg.norm(f(U))
            while np.linalg.norm(f(U + alpha * delta)) > norme :
                alpha *= rho
                if alpha < 1e-10:  # Éviter les boucles infinies
                    break
                U_new = U + alpha * delta
        U = U_new
        f_u.append(np.linalg.norm(f(U)))
    return U, f_u


def test_Newton_Raphson():
    # Define the function f and its derivative J for one dimension
    def f(U):
        x = U[0]
        return np.array([x**2 - 4])

    def J(U):
        x = U[0]
        return np.array([[2*x]])

    # Initial guess
    U0 = np.array([1.0])
    # Number of iterations
    N = 10
    # Tolerance
    epsilon = 1e-6

    # Expected solution
    expected_solution = np.array([2.0])

    # Run the Newton-Raphson method
    solution, _ = Newton_Raphson(f, J, U0, N, epsilon)
    print(f"solution {solution[0]} / expected solution {expected_solution}")

    # Check if the solution is close to the expected solution
    assert np.allclose(solution, expected_solution, atol=epsilon), f"Expected {expected_solution}, but got {solution}"

def test_Newton_Raphson_2D():
    # Define the function f and its Jacobian J for two dimensions
    def f(U):
        x, y = U
        return np.array([x**2 + y**2 - 4, x - y])

    def J(U):
        x, y = U
        return np.array([[2*x, 2*y], [1, -1]])

    # Initial guess
    U0 = np.array([1.0, 1.0])
    # Number of iterations
    N = 10
    # Tolerance
    epsilon = 1e-6

    # Expected solution (one possible solution)
    expected_solution = np.array([np.sqrt(2), np.sqrt(2)])

    # Run the Newton-Raphson method
    solution, _ = Newton_Raphson(f, J, U0, N, epsilon)
    print(f"solution {solution} / expected solution {expected_solution}")

    # Check if the solution is close to the expected solution
    assert np.allclose(solution[0], expected_solution, atol=epsilon), f"Expected {expected_solution}, but got {solution[0]}"

# Run the test
if __name__ == "__main__":
    test_Newton_Raphson()
    test_Newton_Raphson_2D()
