import numpy as np
import matplotlib.pyplot as plt
from matrice import *
from chaleur import *
import time

def test_conjgrad(precisions, num_tests=20):
    sizes = [5 * x for x in range(1, 21)]
    times = {precision: [] for precision in precisions}
    precisions_results = {precision: [] for precision in precisions}
    for size in sizes:
        for precision in precisions:
            time_results = []
            precision_results = []
            for _ in range(num_tests):
                A = mat_sym_pos_def(size, size)
                b = np.random.rand(size)
                x0 = np.zeros(size)
                begin = time.time()
                x = conjgrad(A, b, x0, precision)
                end = time.time() - begin
                time_results.append(end)
                precision_result = np.linalg.norm(np.dot(A, x) - b)
                precision_results.append(precision_result)
                assert np.allclose(np.dot(A, x), b, atol=precision), f"Test failed for size {size} with precision {precision}: The solution does not satisfy A*x = b"
            times[precision].append(np.mean(time_results))
            precisions_results[precision].append(np.mean(precision_results))
    return sizes, times, precisions_results

def test_conjgrad_precondition(precisions, num_tests=20, use_inverse=False, use_cholesky=False):
    sizes = [5 * x for x in range(1, 21)]
    times = {precision: [] for precision in precisions}
    precisions_results = {precision: [] for precision in precisions}
    for size in sizes:
        for precision in precisions:
            time_results = []
            precision_results = []
            for _ in range(num_tests):
                A = mat_sym_pos_def(size, size)
                b = np.random.rand(size)
                x0 = np.zeros(size)
                begin = time.time()
                if use_inverse:
                    M = np.linalg.inv(A)
                elif use_cholesky:
                    M = cholesky_preconditioner(A)
                else:
                    M = np.diag(1 / np.diag(A))
                x = conjgrad_precondition(A, b, x0, M, precision)
                end = time.time() - begin
                time_results.append(end)
                precision_result = np.linalg.norm(np.dot(A, x) - b)
                precision_results.append(precision_result)
                assert np.allclose(np.dot(A, x), b, atol=precision), f"Test failed for size {size} with precision {precision}: The solution does not satisfy A*x = b"
            times[precision].append(np.mean(time_results))
            precisions_results[precision].append(np.mean(precision_results))
    return sizes, times, precisions_results

def plot_times(sizes, times_conjgrad, times_conjgrad_precondition, times_conjgrad_precondition_inv, times_conjgrad_precondition_cholesky, precision):
    plt.plot(sizes, times_conjgrad[precision], label='Conjugate Gradient')
    plt.plot(sizes, times_conjgrad_precondition[precision], label='Conjugate Gradient with Jacobi Preconditioning')
    plt.plot(sizes, times_conjgrad_precondition_inv[precision], label='Conjugate Gradient with A^-1 Preconditioning')
    plt.plot(sizes, times_conjgrad_precondition_cholesky[precision], label='Conjugate Gradient with Cholesky Preconditioning')
    plt.xlabel('Size of the array')
    plt.ylabel('Time (seconds)')
    plt.title(f'Comparison of Conjugate Gradient Methods (precision={precision})')
    plt.legend()
    plt.show()

def plot_precision(sizes, precisions_conjgrad, precisions_conjgrad_precondition, precisions_conjgrad_precondition_inv, precisions_conjgrad_precondition_cholesky, precision):
    plt.plot(sizes, precisions_conjgrad[precision], label='Conjugate Gradient')
    plt.plot(sizes, precisions_conjgrad_precondition[precision], label='Conjugate Gradient with Jacobi Preconditioning')
    plt.plot(sizes, precisions_conjgrad_precondition_inv[precision], label='Conjugate Gradient with A^-1 Preconditioning')
    plt.plot(sizes, precisions_conjgrad_precondition_cholesky[precision], label='Conjugate Gradient with Cholesky Preconditioning')
    plt.xlabel('Size of the array')
    plt.ylabel('Precision (norm of residual)')
    plt.title(f'Precision Comparison of Conjugate Gradient Methods (precision={precision})')
    plt.legend()
    plt.show()
    
def test_res_eq_chaleur():
    N = 20
    F_constant = np.ones(N**2)  # Source de chaleur constante
    F_random = np.random.rand(N**2)  # Source de chaleur aléatoire
    F_center_hot = np.zeros(N**2)  # Source de chaleur au centre
    F_center_hot[(N**2 +N)// 2 ] = 20
    
    print("Test avec source de chaleur constante et méthode itérative:")
    res_eq_chaleur_iter(F_constant) 
    print("Test avec source de chaleur constante et méthode directe:")
    res_eq_chaleur_cholesky(F_constant)
    print("Test avec source de chaleur aléatoire et méthode itérative:")
    res_eq_chaleur_iter(F_random)
    print("Test avec source de chaleur aléatoire et méthode directe:")
    res_eq_chaleur_cholesky(F_random)
    print("Test avec source de chaleur au centre et méthode itérative:")
    res_eq_chaleur_iter(F_center_hot)
    print("Test avec source de chaleur au centre et méthode directe:")
    res_eq_chaleur_cholesky(F_center_hot)


if __name__ == "__main__":
    """
    precisions = [1e-5, 1e-8, 1e-10]
    constant_precision = 1e-8  # Choose a constant precision for the time plot
    sizes_conjgrad, times_conjgrad, precisions_conjgrad = test_conjgrad(precisions)
    sizes_conjgrad_precondition, times_conjgrad_precondition, precisions_conjgrad_precondition = test_conjgrad_precondition(precisions)
    sizes_conjgrad_precondition_inv, times_conjgrad_precondition_inv, precisions_conjgrad_precondition_inv = test_conjgrad_precondition(precisions, use_inverse=True)
    sizes_conjgrad_precondition_cholesky, times_conjgrad_precondition_cholesky, precisions_conjgrad_precondition_cholesky = test_conjgrad_precondition(precisions, use_cholesky=True)
    plot_times(sizes_conjgrad, times_conjgrad, times_conjgrad_precondition, times_conjgrad_precondition_inv, times_conjgrad_precondition_cholesky, constant_precision)
    plot_precision(sizes_conjgrad, precisions_conjgrad, precisions_conjgrad_precondition, precisions_conjgrad_precondition_inv, precisions_conjgrad_precondition_cholesky, constant_precision)
    """
    test_res_eq_chaleur()
    print("All tests passed.")
