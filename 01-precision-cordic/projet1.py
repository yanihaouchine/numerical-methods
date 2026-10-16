import matplotlib.pyplot as plt
import numpy as np
import math

# Question 1

def reduce(x, p):
    """
    Réduit un nombre x à p chiffres significatifs.
    """
    i = 1 
    negative = 0
    if x < 0:
        negative = 1
        x = abs(x)
    if x == 0.0: 
        return x
    if x > 1:
        while x >= 1:
            x /= 10
            i *= 10
    elif x < 0.1:
        while x < 0.1:
            x *= 10 
            i /= 10
    x = round(x, p)
    x *= i
    if negative:
        return -x
    else:
        return x

def test_reduce():
    """
    Teste la fonction reduce avec des valeurs et des précisions données.
    """
    val = [3.141592658, 10507.1823, 0.0001857563]
    result4 = [3.142, 10510, 0.0001858]
    result6 = [3.14159, 10507.2, 0.000185756]
    for i in range(3):
        print("Test pour ", val[i], "sur 4, attendu: ", result4[i], "retourné: ", reduce(val[i], 4))
        print("Test pour ", val[i], "sur 6, attendu: ", result6[i], "retourné: ", reduce(val[i], 6))

test_reduce()

# Question 2 

def add(x, y, p):
    """
    Additionne deux nombres x et y en les réduisant à p chiffres significatifs.
    """
    return reduce(x + y, p)

def mul(x, y, p):
    """
    Multiplie deux nombres x et y en les réduisant à p chiffres significatifs.
    """
    return reduce(x * y, p)

# Question 3 

def delta_sum(x, y, p):
    """
    Calcule l'erreur relative de l'addition de x et y en représentation décimale réduite.
    """
    x_r = reduce(x, p)
    y_r = reduce(y, p)
    sum_reduced = reduce(x_r + y_r, p)
    sum_real = x + y
    return abs(sum_real - sum_reduced) / abs(sum_real)

# Question 4

def delta_mul(x, y, p):
    """
    Calcule l'erreur relative de la multiplication de x et y en représentation décimale réduite.
    """
    x_r = reduce(x, p)
    y_r = reduce(y, p)
    return abs((x * y) - mul(x_r, y_r, p)) / abs(x * y)

# Fixer x et p
x_fixed = 5
# Générer des valeurs pour y
y_values = np.linspace(-10, 10, 200)  # 200 points entre -10 et 10

# Différentes valeurs de précision p
precisions = [4, 5, 6]

# Tracer les courbes pour la somme
plt.figure(figsize=(10, 6))
for p in precisions:
    z_values_s = [delta_sum(x_fixed, y, p) for y in y_values]
    plt.plot(y_values, z_values_s, label=f"$p={p}$")

# Ajouter des détails au graphique pour la somme
plt.title("Erreur relative de la somme en fonction de $y$", fontsize=14)
plt.xlabel("$y$", fontsize=12)
plt.ylabel("$\\delta_{s}(x, y, p)$", fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)

colors = ["blue", "orange", "green"]
for i in range(3):
    plt.hlines(10**(-precisions[i]+1)/2, -10, 10, colors[i])

# Afficher le graphique pour la somme
plt.show()

# Tracer les courbes pour la multiplication
plt.figure(figsize=(10, 6))
for p in precisions:
    z_values_m = [delta_mul(x_fixed, y, p) for y in y_values]
    plt.plot(y_values, z_values_m, label=f"$p={p}$")

# Ajouter des détails au graphique pour la multiplication
plt.title("Erreur relative de la multiplication en fonction de $y$", fontsize=14)
plt.xlabel("$y$", fontsize=12)
plt.ylabel("$\\delta_{p}(x, y, p)$", fontsize=12)
plt.legend(fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)

colors = ["blue", "orange", "green"]
for i in range(3):
    plt.hlines(10**(-precisions[i]+1)/2, -10, 10, colors[i])

# Afficher le graphique pour la multiplication
plt.show()

# PARTIE 2
# Question 1
import math

def log2_approx(p):
    """
    Approxime log(2) à p décimales en utilisant la série de Taylor.
    """
    log2 = 0
    n = 1
    term = 1
    while abs(term) > 10**(-p):
        if (n % 2) == 0:
            term = -1 / n
        if (n % 2) == 1:
            term = 1 / n
        log2 += term
        n += 1
    
    true_log2 = math.log(2)
    relative_error = abs((true_log2 - log2) / true_log2)
    
    return round(log2, p), relative_error

# Exemple
p = 4 
approx_log2, error = log2_approx(p)
print(f"Valeur approchée de log(2) sur {p} décimales : {approx_log2}")
print(f"Erreur relative : {error}")

# Question 2

# Pré-calcul des valeurs de log(1 + 10^(-k)) et atan(10^(-k))
L = [math.log(1 + 10**(-k)) for k in range(6)]
A = [math.atan(10**(-k)) for k in range(4)]

def ln(x):
    """
    Calcule ln(x) en utilisant une méthode basée sur les valeurs pré-calculées de log(1 + 10^(-k)).
    """
    k = 0
    y = 0
    p = 1
    while (k < 6):
        while (x >= p + p * 10**(-k)):
            y = y + L[k]
            p = p + p * 10**(-k)
        k += 1
    return y + (x / p - 1)

def exp(x):
    """
    Calcule exp(x) en utilisant une méthode basée sur les valeurs pré-calculées de log(1 + 10^(-k)).
    """
    k = 0
    y = 1
    while (k < 6):
        while (x >= L[k]):
            x = x - L[k]
            y = y + y * 10**(-k)
        k += 1
    return y + y * x

def arctan(x):
    """
    Calcule arctan(x) en utilisant une méthode basée sur les valeurs pré-calculées de atan(10^(-k)).
    """
    k = 0
    y = 1
    r = 0
    while (k < 4):
        while (x < y * 10**(-k)):
            k += 1
        xp = x - y * 10**(-k)
        y = y + x * 10**(-k)
        x = xp
        r = r + math.atan(10**(-k))
    return r + (x / y)

def tan(x):
    """
    Calcule tan(x) en utilisant une méthode basée sur les valeurs pré-calculées de atan(10^(-k)).
    """
    k = 0
    n = 0
    d = 1
    for k in range(4):
        while (x >= A[k]):
            x = x - A[k]
            np = n + d * 10**(-k)
            d = d - n * 10**(-k)
            n = np
    return (n + x * d) / (d - x * n)

# Tests des fonctions
print("exp(1) =", exp(1))  # Devrait être proche de 2.7183
print("ln(2) =", ln(2))    # Devrait être proche de 0.69314
print("tan(1) =", tan(1))  # Devrait être proche de 1.5574
print("arctan(1) =", arctan(1))  # Devrait être proche de 0.7854