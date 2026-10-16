import numpy as np

def HH_vector(x):
    """Construit un vecteur de Householder
    
    Args:
        x (ndarray): Vecteur de départ.

    Returns:
        ndarray: Vecteur de Householder v.
    """
    v = x.copy()
    v[0] += np.sign(x[0]) * np.linalg.norm(x)
    v /= np.linalg.norm(v)
    return v

def householder_matrix(U, V):
    """Construit la matrice de Householder qui transforme U en V.

    Args:
        U (ndarray): Vecteur de départ.
        V (ndarray): Vecteur d'arrivée.

    Returns:
        ndarray: Matrice de Householder H.
    """
    N = U - V
    N = N / np.linalg.norm(N)  # Normalisation pour rendre N unitaire
    I = np.eye(len(U))  # Matrice identité
    H = I - 2 * np.outer(N, N.T)  # Construction de H
    return H

def HH_matrix(N):
    """Construit la matrice de Householder à partir d'un vecteur unitaire N.

    Args:
        N (ndarray): Vecteur unitaire.

    Returns:
        ndarray: Matrice de Householder H.
    """
    I = np.eye(len(N))
    H = I - 2 * np.outer(N, N.T)
    return H

def apply_householder_to_vector(H, U):
    """Applique la matrice de Householder H à un vecteur U.

    Args:
        H (ndarray): Matrice de Householder.
        U (ndarray): Vecteur à transformer.

    Returns:
        ndarray: Vecteur transformé.
    """
    return H @ U

def apply_householder_to_matrix(H, A):
    """Applique la matrice de Householder H à une matrice A (colonne par colonne).

    Args:
        H (ndarray): Matrice de Householder.
        A (ndarray): Matrice à transformer.

    Returns:
        ndarray: Matrice transformée.
    """
    return H @ A


if __name__ == '__main__':
    # Définition des vecteurs U et V
    U = np.array([[3], [4], [0]])
    V = np.array([[0], [0], [5]])

    # Construction de la matrice de Householder
    H = householder_matrix(U, V)

    # Application de H à U
    np.set_printoptions(suppress=True)
    U_transformed = apply_householder_to_vector(H, U)
    print("Matrice de Householder H :\n", H)
    print("\nH * U =\n", U_transformed)
    print("\nV attendu :\n", V)

    # Vérification de la transformation
    assert np.allclose(U_transformed, V), "Erreur : H(U) n'est pas égal à V"

    # Construction de la matrice de Householder à partir de N
    N = U - V
    N = N / np.linalg.norm(N)
    H = HH_matrix(N)
    print("\nMatrice de Householder H :\n", H)

    # Application de H à une matrice
    A = np.array([[1, 2], [3, 4], [5, 6]])
    A_transformed = apply_householder_to_matrix(H, A)
    print("\nH * A =\n", A_transformed)
    print("\nA attendu :\n", A)
    print("\nVérification : H * A =\n", H @ A)
    assert np.allclose(A_transformed, H @ A), "Erreur : H(A) n'est pas égal à A"
