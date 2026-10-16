import numpy as np
from matrice import *


def cholesky_decomposition(A):
    n = A.shape[0]
    T = np.zeros_like(A)
    for i in range(n):
        for j in range(i + 1):
            sum_k = sum(T[i, k] * T[j, k] for k in range(j))            
            if i == j:
                T[i, j] = np.sqrt(A[i, i] - sum_k)
            else:
                T[i, j] = (A[i, j] - sum_k) / T[j, j]
    
    return T

# Tester l'algorithme
A = np.array([[4, 12, -16], 
              [12, 37, -43], 
              [-16, -43, 98]], dtype=float)

T = cholesky_decomposition(A)


  
#print("Matrice T obtenue par décomposition de Cholesky:")
#print(T)

def transposé(T):
    n=T.shape[0]
    M=np.zeros_like(T)
    for i in range(n):
        for j in range(n):
            M[i][j]=T[j][i]
    return M

#Tester la factorisation complète :
print("La transposé de T est: ")
T_transposé=transposé(T)
print(T_transposé)

print("La matrice A est:")
print(A)

print("Le résultat de T multiplié par la transposé de T est :")
résultat=np.dot(T,T_transposé)
print(résultat)






#Factorisation incomplète

def factorisation_incomplete(A):
    n = A.shape[0]
    T = np.zeros_like(A)  

    for i in range(n):
        for j in range(i + 1):
            if A[i, j] != 0:  
                sum_k = sum(T[i, k] * T[j, k] for k in range(j))
                if i == j:
                    T[i, j] = np.sqrt(A[i, i] - sum_k)
                else:
                    T[i, j] = (A[i, j] - sum_k) / T[j, j]
    
    return T



#tester la fctorisation incomplète

print("La matrice symétrique creuse est: ")
M=mat_sym_pos_def(4,9)
print(M)

T2=factorisation_incomplete(M)
print("Matrice T obtenue par décomposition incomplète de Cholesky:")
print(T)

T2_transposé=transposé(T2)
résultat=np.dot(T2,T2_transposé)
print("Le résultat de T multiplié par la transposé de T est :")
print(résultat)


#premier cas:
T1=np.zeros_like(M)
T1=cholesky_decomposition(M)
T1_inv=np.linalg.inv(T1)
transposé_T1_inv=transposé(T1_inv)
résultat=np.dot(transposé_T1_inv,T1_inv)
COND_1=np.linalg.cond(résultat)
print("Préconditionnement avec factorisation complète:")
print(COND_1)


#Deuxième cas
T2_inv=np.linalg.inv(T2)
transposé_T2_inv=transposé(T2_inv)
résultat=np.dot(transposé_T2_inv,T2_inv)
COND_2=np.linalg.cond(résultat)
print("Préconditionnement avec factorisation incomplète:")
print(COND_2)


