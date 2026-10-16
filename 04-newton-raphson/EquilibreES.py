import numpy as np
import matplotlib.pyplot as plt
import NewtonRaphson as NR
import numpy.polynomial.legendre as leg

#Calcul de l'energie E d'un système de N charges (première formule du sujet)
def energie(x):
    N = len(x)
    E = 0
    for i in range(N):
        xi = x[i]
        E += np.log(abs(xi + 1)) + np.log(abs(xi - 1))
        for j in range(N):
            if i != j:
                E += 0.5 * np.log(abs(xi - x[j]))
    return E

#Calcul du gradient de l'energie (deuxième formule du sujet)
def gradient(x):
    N = len(x)
    grad = np.zeros(N)
    for i in range(N):
        xi = x[i]
        gE = 0
        gE += 1 / (xi + 1) + 1 / (xi - 1)
        for j in range(N):
            if i != j:
                gE += 1 / (xi - x[j])
        grad[i] = gE
    return grad

#Calcul de la jacobienne du gradient (question 1 du sujet)
def jacobienne_du_gradient(x):
    N = len(x)
    J = np.zeros((N, N))
    for i in range(N):
        xi = x[i]
        for j in range(N):
            if i != j:
                J[i, j] = 1 / (xi - x[j])**2
            else:
                J[i, i] = -1 / (xi + 1)**2 - 1 / (xi - 1)**2
                for k in range(N):
                    if k != i:
                        J[i, i] -= 1 / (xi - x[k])**2
    return J

def solve_with_NR(N,iter):
    Uc = np.linspace(-0.9,0.9,N)
    print("Vecteur U initial représentant les N charges :")
    print(Uc)
    U, diff = NR.Newton_Raphson(gradient,jacobienne_du_gradient,Uc,iter,1e-6)
    _, diff_no_bt= NR.Newton_Raphson(gradient,jacobienne_du_gradient,Uc,iter,1e-6, backtracking=False)
    return U, diff, diff_no_bt

def plot_points(points) :
    derive_legendre = leg.Legendre.basis(len(points)+1).deriv()
    x = np.linspace(-1, 1, 1000)
    y = derive_legendre(x) 
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, label=f"Dérivée des polynomes de Legendre de degré N+1", color='r')
    plt.scatter(points, np.zeros_like(points), color='green', label='Charges')
    plt.axvline(-1, color='b', linestyle='--')
    plt.axvline(1, color='b', linestyle='--')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.legend()
    plt.xlabel("Position des charges")
    plt.grid(True)
    plt.show()

def est_definie_positive(A):
    if not np.allclose(A, A.T):
        return False
    
    valeurs_propres = np.linalg.eigvals(A)

    return np.all(valeurs_propres > 0)

def est_definie_negative(A):
    if not np.allclose(A, A.T):
        return False
    
    valeurs_propres = np.linalg.eigvals(A)

    return np.all(valeurs_propres < 0)

#===========================================================================
#================================== TESTS ==================================
#===========================================================================

def test_energie():
    x_test = np.array([-0.5, 0, 0.5]) 
    E = energie(x_test) 
    print("Energie pour les positions [-0.5, 0, 0.5]: \n", E)

def test_gradient():
    x_test = np.array([-0.5, 0, 0.5])
    grad = gradient(x_test) 
    print("Gradient pour les positions [-0.5, 0, 0.5]: \n", grad)

def test_jacobienne_du_gradient():
    x_test = np.array([-0.5, 0, 0.5]) 
    J = jacobienne_du_gradient(x_test) 
    print("Jacobienne pour les positions [-0.5, 0, 0.5]: \n")
    print(J)

def test():
    print("----- Test de la fonction énergie -----")
    test_energie()
    
    print("\n----- Test de la fonction gradient -----")
    test_gradient()
    
    print("\n----- Test de la fonction jacobienne du gradient -----")
    test_jacobienne_du_gradient()

if __name__ == "__main__":
    #test()

    points, diff, diff_no_bt = solve_with_NR(5,100)
    print("Vecteur associé aux positions des N charges pour lesquelles l'énergie totale atteint une valeur extrême :")
    print(points)
    plot_points(points)

    # Plot the convergence speed
    plt.figure(figsize=(8, 6))
    plt.semilogy(range(len(diff_no_bt)), diff_no_bt, 'r-o', label='Without Backtracking')
    plt.semilogy(range(len(diff)), diff, 'b-o', label='With Backtracking')
    plt.xlabel('Iteration')
    plt.ylabel('f(U) (log scale)')
    plt.title('Convergence Speed for Equilibrium with and without Backtracking')
    plt.legend()
    plt.grid()
    plt.show()

    print("Valeurs propres de la jacobienne associée au gradient de l'énergie de ces N charges :")
    print(np.linalg.eigvals(jacobienne_du_gradient(points)))
    print("La jacobienne du gradient d'énergie est définie positive ?", est_definie_positive(jacobienne_du_gradient(points)))
    print("La jacobienne du gradient d'énergie est définie négative ?", est_definie_negative(jacobienne_du_gradient(points)))
    # La jacobienne du gradient de l'énergie, qui est la hessienne de l'énergie, est définie négative. Donc, l'énergie associée aux N charges est maximal pour la position des N charges trouvée précédemment (positions associées au vecteur 'points').

    
#===========================================================================
