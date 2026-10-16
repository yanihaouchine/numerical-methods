from matrice import *
import numpy as np
import matplotlib.pyplot as plt
from cholesky import *
from scipy.linalg import solve_triangular



def create_laplacien_matrix(N):
    # Créer la matrice A représentant le Laplacien
    # Créer un bloc de taille N*N avec -4 sur la diagonale
    # 1 sur les diagonales voisines, 1 sur les diagonales de distance N+1 avec la diagonale centrale
    A=np.zeros((N**2,N**2))
    for i in range(N**2):
        A[i,i]=4
        if i - 1 >= 0:  # Diagonale au-dessus de la diagonale principale
            A[i, i - 1] =- 1
        if i + 1 < N**2:  # Diagonale en-dessous de la diagonale principale
            A[i, i + 1] = -1
        if i - N-1 >= 0:  # Diagonale au-dessus à une distance N+1 de la diagonale
            A[i, i - N - 1] = -1
        if i + N - 1 <= N**2:  # Diagonale en-dessous à une distance N+1 de la diagonale
            A[i - N - 1, i] = -1
        if i%(N-1)==0:
            A[i,i-1]=0
            A[i-1,i]=0           
    print(A)
    return A
        
    
    

def heatmap(T):
    # Créé une heatmap pour représenter plus visuellement l'évolution de la chaleur
    # T : vecteur de taille N**2
    N = int(np.sqrt(len(T)))  # Calculer la taille de la grille
    if N * N != len(T):
        raise ValueError("La taille du vecteur T doit être un carré parfait.")
    
    grid = T.reshape((N, N))  # Reshaper le vecteur en une matrice 2D
    
    plt.imshow(grid, cmap='hot', interpolation='nearest', origin='lower')
    plt.colorbar(label='Temperature')
    plt.title('Heatmap de la distribution de chaleur')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
    
def res_eq_chaleur(F):
    N = int(np.sqrt(len(F)))
    A = create_laplacien_matrix(N)
    
    # Initialisation de la température initiale à zéro
    Temperature = np.zeros(N**2)
    
    # Résout l'équation de la chaleur avec la méthode du gradient conjugué
    T = conjgrad(A, F, Temperature)
    """
    # Avec cholesky
    L=cholesky_decomposition(A)
    L=-L
    y=solve_triangular(L,F)
    T=solve_triangular(np.transpose(L),y,upper=True)
    """
    print("Solution :")
    heatmap(T)
    
    # Vérifie que la solution est bonne 
    # Utilisation de la norme infinie pour s'assurer que chaque composante est assez précise.
    reste = np.linalg.norm(np.dot(A, T) - F, ord=np.inf)
    print(f"Imprécision: {reste}")
    #assert reste < 1e-5, "The solution does not meet the required precision"

