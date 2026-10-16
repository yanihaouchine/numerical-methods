import numpy as np
 
def mat_sym_pos_def(n, length):
    # Initialise une matrice symétrique définie positive composée de n (ou n+1 selon la parité)
    # éléments extra-diagonaux non nuls
    # n : int -> nombre d'éléments extra-diagonaux non nuls
    #length : int -> taille de la matrice carrée créée
    A = np.zeros((length, length))
    nb_extra_dag = 0
    while nb_extra_dag < n:
        i, j = np.random.randint(0, length, 2)
        while i == j:
            i, j = np.random.randint(0, length, 2)
        if A[i, j] == 0:
            A[i, j] = abs(np.random.rand())
            A[j, i] = A[i, j]
            nb_extra_dag += 2
    for i in range(length):
        sum = 0
        for j in range(length):
            sum += A[i, j]
        A[i, i] = abs(sum) + 1
    return A

def is_pos_def(A):
    # Vérifie que A est symétrique définie positive
    if not np.allclose(A, A.T):
        return False
    return np.all(np.linalg.eigvals(A) > 0)

def conjgrad(A, b, x, precision=1e-10):
    # Résout le système linéaire Ax=b par la méthode du gradient conjugué
    # A : Matrice de taille len(x)*len(x)
    # b : vecteur de taille len(x)
    # x : un vecteur ( on l'initialise arbitrairement à 0 )
    # precision : modifiale selon l'éxigence attendue sur la précision de la résolution
    # retourne x modifié pour correspondre à la solution
    r = b - np.dot(A, x)
    p = r
    rsold = np.dot(r.T, r)

    for i in range(len(b)):
        Ap = np.dot(A, p)
        alpha = rsold / np.dot(p.T, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rsnew = np.dot(r.T, r)
        if np.sqrt(rsnew) < precision:
            break
        p = r + (rsnew / rsold) * p
        rsold = rsnew

    return x

def conjgrad_precondition(A, b, x, M, precision=1e-10):
    # Résout le système linéaire Ax=b par la méthode du gradient conjugué couplé au préconditionnement de A
    # A : Matrice de taille len(x)*len(x)
    # b : vecteur de taille len(x)
    # x : un vecteur ( on l'initialise arbitrairement à 0 )
    # M : Matrice de préconditionnement de A
    # precision : modifiale selon l'éxigence attendue sur la précision de la résolution
    # retourne x modifié pour correspondre à la solution
    r = b - np.dot(A, x)
    z = np.dot(M, r)
    p = z
    rsold = np.dot(r.T, z)

    for i in range(len(b)):
        Ap = np.dot(A, p)
        alpha = rsold / np.dot(p.T, Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        z = np.dot(M, r)
        rsnew = np.dot(r.T, z)
        if np.sqrt(rsnew) < precision:
            break
        p = z + (rsnew / rsold) * p
        rsold = rsnew

    return x

"""
def cholesky_preconditioner(A):
    # Temporaire
    L = np.linalg.cholesky(A)
    L_inv = np.linalg.inv(L)
    return np.dot(L_inv.T, L_inv)
"""