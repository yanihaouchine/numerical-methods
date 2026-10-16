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
    return reduce(reduce(x,p) + reduce(y, p), p)  #modif

def mul(x, y, p):
    """
    Multiplie deux nombres x et y en les réduisant à p chiffres significatifs.
    """
    return reduce(reduce(x,p) * reduce(y, p), p)  #modif

# Question 3 

def delta_sum(x, y, p):
    """
    Calcule l'erreur relative de l'addition de x et y en représentation décimale réduite.
    """
    sum_reduced = add(x,y,p)
    sum_real = x + y
    return abs(sum_real - sum_reduced) / abs(sum_real) #modif

# Question 4

def delta_mul(x, y, p):
    """
    Calcule l'erreur relative de la multiplication de x et y en représentation décimale réduite.
    """
    prod_reduced = mul(x,y,p)
    prod_real = x*y
    return abs(prod_real - prod_reduced) / abs(prod_real) #modif

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


def log2_approx(p):   #modif
    """
    Approxime log(2) à p décimales en utilisant la série de Taylor.
    """
    log2 = 1
    new_log2 = 0
    n = 1
    term = 1
    stop = 10 **(-p)
    while ((log2) != new_log2)and(abs(1/n) > stop):
        log2 = new_log2
        if (n % 2) == 0:
            term = reduce(-1 / n, p)
        if (n % 2) == 1:
            term = reduce(1 / n, p)
        new_log2 = add(new_log2,term,p)
        n += 1
    
    true_log2 = math.log(2)
    relative_error = abs((true_log2 - log2) / true_log2)
    
    return log2, relative_error

# Exemple
p = 4 
approx_log2, error = log2_approx(p)
print(f"Valeur approchée de log(2) sur {p} décimales : {approx_log2}")
print(f"Erreur relative : {error}")

# Question 2

# Pré-calcul des valeurs de log(1 + 10^(-k)) et atan(10^(-k))
L = [math.log(1 + 10**(-k)) for k in range(6)]
A = [math.atan(10**(-k)) for k in range(5)]

def ln_cordic(x):
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

def exp_cordic(x):
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

def atan_cordic(x):  
    k = 0
    y = 1
    r = 0
    
    while k <= 4:
        while x < y * 10**-k:
            k += 1
            if k > 4:
                break  
        
        if k > 4:
            break
        
        xp = x - y * 10**-k
        y = y + x * 10**-k
        x = xp
        r += A[k]
    
    return r + x / y

def tan_cordic(x):
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
print("exp(1) =", exp_cordic(1))  # Devrait être proche de 2.7183
print("ln(2) =", ln_cordic(2))    # Devrait être proche de 0.69314
print("tan(1) =", tan_cordic(1))  # Devrait être proche de 1.5574
print("arctan(1) =", atan_cordic(1))  # Devrait être proche de 0.7854
# Plage de valeurs pour comparer les fonctions
x_values = np.linspace(0.1, 2, 100)  # Éviter x=0 pour log

# Calcul des valeurs
ln_values = [ln_cordic(x) for x in x_values]
exp_values = [exp_cordic(x) for x in x_values]
arctan_values = [atan_cordic(x) for x in np.linspace(-1, 1, 100)]
tan_values = [tan_cordic(x) for x in np.linspace(-1, 1, 100)]

# Références des bibliothèques
ln_lib = np.log(x_values)
exp_lib = np.exp(x_values)
arctan_lib = np.arctan(np.linspace(-1, 1, 100))
tan_lib = np.tan(np.linspace(-1, 1, 100))

# Tracé des fonctions
plt.figure(figsize=(10, 6))

plt.subplot(2, 2, 1)
plt.plot(x_values, ln_values, label='ln(x) (custom)', linestyle='dashed')
plt.plot(x_values, ln_lib, label='ln(x) (numpy)')
plt.legend()
plt.title("Comparaison ln(x)")

plt.subplot(2, 2, 2)
plt.plot(x_values, exp_values, label='exp(x) (custom)', linestyle='dashed')
plt.plot(x_values, exp_lib, label='exp(x) (numpy)')
plt.legend()
plt.title("Comparaison exp(x)")

plt.subplot(2, 2, 3)
plt.plot(np.linspace(-1, 1, 100), arctan_values, label='arctan(x) (custom)', linestyle='dashed')
plt.plot(np.linspace(-1, 1, 100), arctan_lib, label='arctan(x) (numpy)')
plt.legend()
plt.title("Comparaison arctan(x)")

plt.subplot(2, 2, 4)
plt.plot(np.linspace(-1, 1, 100), tan_values, label='tan(x) (custom)', linestyle='dashed')
plt.plot(np.linspace(-1, 1, 100), tan_lib, label='tan(x) (numpy)')
plt.legend()
plt.title("Comparaison tan(x)")

plt.tight_layout()
plt.show()

def erreur_relative(approx, exact):
    return np.abs((approx - exact) / exact)

# Définition des valeurs de x
x_ln_exp = np.linspace(0.1, 2, 100)  # Intervalle pour ln et exp
x_trig = np.linspace(-0.9, 0.9, 100)  # Intervalle pour arctan et tan

# Calcul des erreurs relatives
err_ln = erreur_relative(np.array([ln_cordic(x) for x in x_ln_exp]), np.log(x_ln_exp))
err_exp = erreur_relative(np.array([exp_cordic(x) for x in x_ln_exp]), np.exp(x_ln_exp))
err_arctan = erreur_relative(np.array([atan_cordic(x) for x in x_trig]), np.arctan(x_trig))
err_tan = erreur_relative(np.array([tan_cordic(x) for x in x_trig]), np.tan(x_trig))

# Tracé des erreurs relatives sur la même échelle
y_max = max(np.max(err_ln), np.max(err_exp), np.max(err_arctan), np.max(err_tan)) * 100

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Graphique 1 : Erreur relative de ln(x)
axes[0, 0].plot(x_ln_exp, err_ln * 100, label="Erreur rel. ln(x)", linestyle="--")
axes[0, 0].set_ylim(0, y_max)
axes[0, 0].set_title("Erreur relative ln(x) (%)")
axes[0, 0].legend()

# Graphique 2 : Erreur relative de exp(x)
axes[0, 1].plot(x_ln_exp, err_exp * 100, label="Erreur rel. exp(x)", linestyle="--")
axes[0, 1].set_ylim(0, y_max)
axes[0, 1].set_title("Erreur relative exp(x) (%)")
axes[0, 1].legend()

# Graphique 3 : Erreur relative de arctan(x)
axes[1, 0].plot(x_trig, err_arctan * 100, label="Erreur rel. arctan(x)", linestyle="--")
axes[1, 0].set_ylim(0, y_max)
axes[1, 0].set_title("Erreur relative arctan(x) (%)")
axes[1, 0].legend()

# Graphique 4 : Erreur relative de tan(x)
axes[1, 1].plot(x_trig, err_tan * 100, label="Erreur rel. tan(x)", linestyle="--")
axes[1, 1].set_ylim(0, y_max)
axes[1, 1].set_title("Erreur relative tan(x) (%)")
axes[1, 1].legend()

plt.tight_layout()
plt.show()

def ln_cordic_approx(x,pr):
    k = 0
    y = 0
    p = 1
    while (k < 6):
        while (x >= p + p * 10**(-k)):
            y = add(y,L[k],pr)
            p = add(p, mul(p, 10**(-k),pr), pr)
        k += 1
    true_log = math.log(x)
    log =  add(y , mul(x, 1/ p,pr) , pr) -1
    relative_error = abs((true_log - log) / true_log)
    return log, relative_error
log, relative = ln_cordic_approx(2,4)
print("Rapport entre l'erreur relative de l'algorithme cordic et de l'erreur de l'algorithme utilsant les séries entiérs pour x = 2 :", error/relative)
