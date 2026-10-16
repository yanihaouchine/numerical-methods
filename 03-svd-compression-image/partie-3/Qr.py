import numpy as np
import matplotlib.pyplot as plt

def double_qr_decomposition(BD, Nmax):
    """
    Implémente l'algorithme Double QR :
      U = Id, V = Id, S = BD
      Pour i = 0,...,Nmax-1 :
          (Q1, R1) = qr(S^T)
          (Q2, R2) = qr(R1^T)
          S = R2
          U = U @ Q2
          V = Q1^T @ V

    Retourne (U, S, V).
    """
    m, n = BD.shape
    
    # Initialisation
    U = np.eye(m)
    V = np.eye(n)
    S = BD.copy()
    
    for i in range(Nmax):
        Q1, R1 = np.linalg.qr(S.T)
        Q2, R2 = np.linalg.qr(R1.T)
        
        S = R2
        U = U @ Q2
        V = Q1.T @ V

    return U, S, V

def off_diagonal_norm(A):
    """
    Calcule la norme (Frobenius) des termes hors-diagonale de A.
    C'est un indicateur de "distance" entre A et sa diagonale.
    """
    # On extrait la diagonale
    diag_part = np.diag(np.diag(A))
    # On calcule la norme Frobenius de la différence A - diag(A)
    return np.linalg.norm(A - diag_part, ord='fro')

def test_and_plot_convergence(BD, Nmax):
    """
    Teste la convergence de l'algorithme Double QR sur la matrice BD,
    et trace l'évolution de :
      - La norme hors-diagonale de S
      - L'erreur d'invariance ||U * S * V - BD||
    à chaque itération.
    """
    m, n = BD.shape
    # Initialisation de U, S, V
    U = np.eye(m)
    V = np.eye(n)
    S = BD.copy()
    
    # Listes pour stocker les mesures à chaque itération
    offdiag_list = []
    invariant_error_list = []
    
    # Mesures initiales (itération 0)
    offdiag_list.append(off_diagonal_norm(S))
    inv_err = np.linalg.norm(U @ S @ V - BD, ord='fro') / np.linalg.norm(BD, ord='fro')
    invariant_error_list.append(inv_err)
    
    # Boucle
    for i in range(Nmax):
        # Etape 1 : decomposition QR sur S^T
        Q1, R1 = np.linalg.qr(S.T)
        # Etape 2 : decomposition QR sur R1^T
        Q2, R2 = np.linalg.qr(R1.T)
        
        # Mise à jour
        S = R2
        U = U @ Q2
        V = Q1.T @ V
        
        # Mesures
        offdiag_list.append(off_diagonal_norm(S))
        inv_err = np.linalg.norm(U @ S @ V - BD, ord='fro') / np.linalg.norm(BD, ord='fro')
        invariant_error_list.append(inv_err)
    
    # Affichage des résultats
    # ========================
    plt.figure(figsize=(8, 5))
    plt.semilogy(offdiag_list, label='Norme hors-diagonale de S')
    plt.semilogy(invariant_error_list, label='Erreur relative ||U S V - BD||')
    plt.xlabel('Itération')
    plt.ylabel('Valeur (échelle log)')
    plt.title('Convergence du Double QR')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return U, S, V, offdiag_list, invariant_error_list

if __name__ == "__main__":
    # Exemple : matrice BD aléatoire 5x5 (carrée pour simplifier le dessin)
    np.random.seed(123)
    BD = np.random.randn(5, 5)
    
    Nmax = 20
    
    U, S, V, offdiag_list, invariant_error_list = test_and_plot_convergence(BD, Nmax)
    
    print("Dernière valeur hors-diagonale de S =", offdiag_list[-1])
    print("Dernière erreur relative ||U S V - BD|| =", invariant_error_list[-1])
    
    # Vérification numérique (exemple)
    reconstruction_error = np.linalg.norm(U @ S @ V - BD, ord='fro')
    print(f"Erreur absolue sur la reconstruction : {reconstruction_error:.2e}")
