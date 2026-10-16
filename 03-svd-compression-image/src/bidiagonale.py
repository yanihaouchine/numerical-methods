import numpy as np
import matplotlib.pyplot as plt
from householder import HH_matrix, HH_vector

norms = []



def mise_sous_forme_bidiagonale(A):
    """
    Bidiagonalisation de la matrice A avec des transformations de Householder.
    
    Args:
        A (ndarray): matrice

    Returns:
        Qleft: Changements de base à gauche
        BD: matrice sous forme bidiagonale
        Qright: Changements de base à droite


    """
    n, m = A.shape
    Qleft = np.eye(n)
    Qright = np.eye(m)
    BD = A.copy()
    
    for i in range(n):
        v1 = HH_vector(BD[i:, i])
        Q1 = np.eye(n)
        Q1[i:, i:] = HH_matrix(v1)
        Qleft = Qleft @ Q1
        BD = Q1 @ BD


        if not np.allclose(Qleft @ BD @ Qright, A): # Vérifier que l’invariant est conservé
            print("Erreur")
            break
        
        if i != m-2:   
            v2 = HH_vector(BD[i, i+1:])
            Q2 = np.eye(m)
            Q2[i+1:, i+1:] = HH_matrix(v2)
            Qright = Q2 @ Qright
            BD = BD @ Q2
        norm_diff = np.linalg.norm((Qleft @ BD @ Qright) - A)
        norms.append(norm_diff)
    
    return Qleft, BD, Qright

if __name__ == '__main__':
    img = np.mean(plt.imread("../img/p3_takeoff_base.png"), axis=2)
    n, m = img.shape
    Qleft, BD, Qright = mise_sous_forme_bidiagonale(img)
    plt.plot(norms)
    plt.title('Convergence de Qleft * BD * Qright vers A')
    plt.xlabel('Itération')
    plt.ylabel('Erreur')
    plt.savefig('../img/erreur_bidiagionale.png')
    plt.show()
