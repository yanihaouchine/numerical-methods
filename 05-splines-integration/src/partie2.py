
def integration_n_left_rectangle(cos,a,b,N):
    h=(b-a)/N
    c=0
    for i in range(N):
        c+=cos(a+i*h)
    return h*c

def integration_n_right_rectangle(cos,a,b,N):
    h=(b-a)/N
    c=0
    for i in range(1,N+1):
        c+=cos(a+i*h)
    return h*c

def integration_n_midle_point(cos,a,b,N):
    h=(b-a)/N
    c=0
    for i in range(N):
        c+=cos(a+ i*h+ h/2)
    return c*h

def integration_n_trapeze(cos,a,b,N):
    h=(b-a)/N   
    c=0
    for i in range(1,N):
        c+=cos(a+i*h)
    return h*((cos(a)+cos(b))/2 +c)

def integration_n_simpson(cos,a,b,N):
    if N % 2 != 0:
        raise ValueError("N doit être pair")
    h=(b-a)/N
    c1=0
    c2=0
    for i in range(1,N,2):
        c1+=cos(a+i*h)
    for j in range(2,N-1,2):
        c2+=cos(a+j*h)
    return (h*(cos(a)+cos(b)+4*c1 + 2*c2))/3


def integration_epsilon(cos,a,b,eps,meth):
    N=1
    if (meth==integration_n_simpson):
        N = 2
    Iapprox_h1=meth(cos,a,b,N)
    while True:
        N=N*2
        Iapprox_h2=meth(cos,a,b,N)
        if (abs(Iapprox_h1- Iapprox_h2)<=eps or N>1000000):
            return Iapprox_h2
        Iapprox_h1=Iapprox_h2



#tester ces methodes pour la fonction cos(x)
import math
a=0
b=math.pi/2
eps=0.0001
N=100

def cos(x):
    return math.cos(x)


resultat_n_left_rectangle=integration_epsilon(cos,a,b,eps,integration_n_left_rectangle)
resultat_n_right_rectangle=integration_epsilon(cos,a,b,eps,integration_n_right_rectangle)
resultat_n_midle_point=integration_epsilon(cos,a,b,eps,integration_n_midle_point)
resultat_n_trapeze=integration_epsilon(cos,a,b,eps,integration_n_trapeze)
resultat_n_simpson=integration_epsilon(cos,a,b,eps,integration_n_simpson)


print("la valeur exacte est :", 1)
print("Integration par la méthode des rectangles à gauche:", resultat_n_left_rectangle,"erreur:",abs(resultat_n_left_rectangle-1))
print("Integration par la méthode des rectangles à droite:", resultat_n_right_rectangle,"erreur:",abs(resultat_n_right_rectangle-1))
print("Integration par la méthode des points milieux:", resultat_n_midle_point,"erreur:",abs(resultat_n_midle_point-1))
print("Integration par la méthode des trapèzes:", resultat_n_trapeze, "erreur:",abs(resultat_n_trapeze-1))
print("Integration par la méthode de Simpson:", resultat_n_simpson, "erreur:",abs(resultat_n_simpson-1))

# >>>>>D'apres les resultats obtenus on peut observer que la méthode de Simpson est la plus précise, suivie de la méthode des trapèzes, puis de la méthode des points milieux, et enfin des méthodes des rectangles à gauche et à droite.

#tester la fonction sin
b=math.pi
def sin(x):
    return math.sin(x)

resultat_n_left_rectangle=integration_epsilon(sin,a,b,eps,integration_n_left_rectangle)
resultat_n_right_rectangle=integration_epsilon(sin,a,b,eps,integration_n_right_rectangle)
resultat_n_midle_point=integration_epsilon(sin,a,b,eps,integration_n_midle_point)
resultat_n_trapeze=integration_epsilon(sin,a,b,eps,integration_n_trapeze)
resultat_n_simpson=integration_epsilon(sin,a,b,eps,integration_n_simpson)

print("la valeur exacte est :", 2)
print("Integration par la méthode des rectangles à gauche:", resultat_n_left_rectangle,"erreur:",abs(resultat_n_left_rectangle-2))
print("Integration par la méthode des rectangles à droite:", resultat_n_right_rectangle,"erreur:",abs(resultat_n_right_rectangle-2))
print("Integration par la méthode des points milieux:", resultat_n_midle_point,"erreur:",abs(resultat_n_midle_point-2))
print("Integration par la méthode des trapèzes:", resultat_n_trapeze, "erreur:",abs(resultat_n_trapeze-2))
print("Integration par la méthode de Simpson:", resultat_n_simpson, "erreur:",abs(resultat_n_simpson-2))


#tester pour la fonction tan
b=math.pi/4
def tan(x):
    return math.tan(x)

resultat_n_left_rectangle=integration_epsilon(tan,a,b,eps,integration_n_left_rectangle)
resultat_n_right_rectangle=integration_epsilon(tan,a,b,eps,integration_n_right_rectangle)
resultat_n_midle_point=integration_epsilon(tan,a,b,eps,integration_n_midle_point)     
resultat_n_trapeze=integration_epsilon(tan,a,b,eps,integration_n_trapeze)
resultat_n_simpson=integration_epsilon(tan,a,b,eps,integration_n_simpson)


print("la valeur exacte est :", 1)
print("Integration par la méthode des rectangles à gauche:", resultat_n_left_rectangle,"erreur:",abs(resultat_n_left_rectangle-1))
print("Integration par la méthode des rectangles à droite:", resultat_n_right_rectangle,"erreur:",abs(resultat_n_right_rectangle-1))
print("Integration par la méthode des points milieux:", resultat_n_midle_point,"erreur:",abs(resultat_n_midle_point-1))
print("Integration par la méthode des trapèzes:", resultat_n_trapeze, "erreur:",abs(resultat_n_trapeze-1))
print("Integration par la méthode de Simpson:", resultat_n_simpson, "erreur:",abs(resultat_n_simpson-1))


#tester pour une fonction polynomiale d'ordre 3
b=1
def pol_4(x):
    return x**4

resultat_n_left_rectangle=integration_epsilon(pol_4,a,b,eps,integration_n_left_rectangle)
resultat_n_right_rectangle=integration_epsilon(pol_4,a,b,eps,integration_n_right_rectangle)
resultat_n_midle_point=integration_epsilon(pol_4,a,b,eps,integration_n_midle_point)
resultat_n_trapeze=integration_epsilon(pol_4,a,b,eps,integration_n_trapeze)
resultat_n_simpson=integration_epsilon(pol_4,a,b,eps,integration_n_simpson)

print("la valeur exacte est :",0,2)
print("Integration par la méthode des rectangles à gauche:", resultat_n_left_rectangle,"erreur:",abs(resultat_n_left_rectangle-0.2))
print("Integration par la méthode des rectangles à droite:", resultat_n_right_rectangle,"erreur:",abs(resultat_n_right_rectangle-0.2))
print("Integration par la méthode des points milieux:", resultat_n_midle_point,"erreur:",abs(resultat_n_midle_point-0.2))
print("Integration par la méthode des trapèzes:", resultat_n_trapeze, "erreur:",abs(resultat_n_trapeze-0.2))
print("Integration par la méthode de Simpson:", resultat_n_simpson, "erreur:",abs(resultat_n_simpson-0.2))
#tester pour une fonction polynomiale d'ordre 4
b=1
def pol_3(x):
    return x**9
resultat_n_left_rectangle=integration_epsilon(pol_3,a,b,eps,integration_n_left_rectangle)
resultat_n_right_rectangle=integration_epsilon(pol_3,a,b,eps,integration_n_right_rectangle)
resultat_n_midle_point=integration_epsilon(pol_3,a,b,eps,integration_n_midle_point)
resultat_n_trapeze=integration_epsilon(pol_3,a,b,eps,integration_n_trapeze)
resultat_n_simpson=integration_epsilon(pol_3,a,b,eps,integration_n_simpson)

print("la valeur exacte est :",0,1)
print("Integration par la méthode des rectangles à gauche:", resultat_n_left_rectangle,"erreur:",abs(resultat_n_left_rectangle-0.1))
print("Integration par la méthode des rectangles à droite:", resultat_n_right_rectangle,"erreur:",abs(resultat_n_right_rectangle-0.1))
print("Integration par la méthode des points milieux:", resultat_n_midle_point,"erreur:",abs(resultat_n_midle_point-0.1))
print("Integration par la méthode des trapèzes:", resultat_n_trapeze, "erreur:",abs(resultat_n_trapeze-0.1))
print("Integration par la méthode de Simpson:", resultat_n_simpson, "erreur:",abs(resultat_n_simpson-0.1))

#partie graphiques*********************************************
import matplotlib.pyplot as plt
a=0
b=math.pi/2
#Camparaison de vitesse de convergence des différentes méthodes d'intégration pour la fonction cos

#on fait varier N qui represente le nombre de subdivisions
N_array=[10,20,50,100,200,300,500,1000]

error_array_left_rectangle=[]
error_array_right_rectangle=[]
error_array_midle_point=[]
error_array_trapeze=[]
error_array_simpson=[]

for i in N_array:
    error_array_left_rectangle.append(abs(integration_n_left_rectangle(cos,a,b,i)-1))
    error_array_right_rectangle.append(abs(integration_n_right_rectangle(cos,a,b,i)-1))
    error_array_midle_point.append(abs(integration_n_midle_point(cos,a,b,i)-1))
    error_array_trapeze.append(abs(integration_n_trapeze(cos,a,b,i)-1))
    error_array_simpson.append(abs(integration_n_simpson(cos,a,b,i)-1))

plt.figure(figsize=(10, 6))
plt.plot(N_array, error_array_left_rectangle, label='Méthode des rectangles à gauche', marker='o')
plt.plot(N_array, error_array_right_rectangle, label='Méthode des rectangles à droite', marker='x')
plt.plot(N_array, error_array_midle_point, label='Méthode des points milieux', marker='d')
plt.plot(N_array, error_array_trapeze, label='Méthode des trapèzes', marker='s')
plt.plot(N_array, error_array_simpson, label='Méthode de Simpson', marker='^')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Erreur absolue de la méthode des rectangles à gauche')
plt.title('Camparaison de vitesse de convergence des méthodes d\'intégration')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()


#Evolution de la valeur de l'intégrale calclulée par la méthode des rectangles à gauche en fonction de N

integrale_array_left_rectangle=[]
integrale_array_right_rectangle=[]
integrale_array_trapeze=[]
integrale_array_simpson=[]
integrale_array_midle_point=[]


for i in N_array:
    integrale_array_left_rectangle.append(integration_n_left_rectangle(cos,a,b,i))
    integrale_array_right_rectangle.append(integration_n_right_rectangle(cos,a,b,i))
    integrale_array_trapeze.append(integration_n_trapeze(cos,a,b,i))
    integrale_array_simpson.append(integration_n_simpson(cos,a,b,i))
    integrale_array_midle_point.append(integration_n_midle_point(cos,a,b,i))

plt.figure(figsize=(10, 6))
plt.plot(N_array, integrale_array_left_rectangle, label='Méthode des rectangles à gauche', marker='o')
plt.axhline(y=1, color='r', label='Valeur exacte de l\'intégrale')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Valeur de l\'intégrale de cos(x) par la méthode des rectangles à gauche')
plt.title('Evolution de la valeur de l\'intégrale par la méthode des rectangles à gauche')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()

#Evolution de la valeur de l'intégrale calclulée par la méthode des rectangles à droite en fonction de N

plt.figure(figsize=(10, 6))
plt.plot(N_array, integrale_array_right_rectangle, label='Méthode des rectangles à droite', marker='x')
plt.axhline(y=1, color='r', label='Valeur exacte de l\'intégrale')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Valeur de l\'intégrale de cos(x) par la méthode des rectangles à droite')
plt.title('Evolution de la valeur de l\'intégrale par la méthode des rectangles à droite')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()
#Evolution de la valeur de l'intégrale calclulée par la méthode des points milieux en fonction de N 
plt.figure(figsize=(10, 6))
plt.plot(N_array, integrale_array_midle_point, label='Méthode des points milieux', marker='d')
plt.axhline(y=1, color='r', label='Valeur exacte de l\'intégrale')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Valeur de l\'intégrale de cos(x) par la méthode des points milieux')
plt.title('Evolution de la valeur de l\'intégrale par la méthode des points milieux')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()
#Evolution de la valeur de l'intégrale calclulée par la méthode des trapèzes en fonction de N  
plt.figure(figsize=(10, 6))
plt.plot(N_array, integrale_array_trapeze, label='Méthode des trapèzes', marker='s')
plt.axhline(y=1, color='r', label='Valeur exacte de l\'intégrale')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Valeur de l\'intégrale de cos(x) par la méthode des trapèzes')
plt.title('Evolution de la valeur de l\'intégrale par la méthode des trapèzes')
plt.grid(True, which='both', linestyle='--', linewidth=0.5) 
plt.legend()
plt.show()

#Evolution de la valeur de l'intégrale calclulée par la méthode de Simpson en fonction de N
plt.figure(figsize=(10, 6))
plt.plot(N_array, integrale_array_simpson, label='Méthode de Simpson', marker='^')
plt.axhline(y=1, color='r', label='Valeur exacte de l\'intégrale')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Valeur de l\'intégrale de cos(x) par la méthode de Simpson')
plt.title('Evolution de la valeur de l\'intégrale par la méthode de Simpson')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()



#Camparaison Evolution de la valeur de l'intégrale calclulée par les différentes méthodes d'intégration en fonction de N 
plt.figure(figsize=(10, 6))
plt.plot(N_array, integrale_array_left_rectangle, label='Méthode des rectangles à gauche', marker='o')
plt.plot(N_array, integrale_array_right_rectangle, label='Méthode des rectangles à droite', marker='x')
plt.plot(N_array, integrale_array_trapeze, label='Méthode des trapèzes', marker='s')
plt.plot(N_array, integrale_array_simpson, label='Méthode de Simpson', marker='^')
plt.plot(N_array, integrale_array_midle_point, label='Méthode des points milieux', marker='d')
plt.axhline(y=1, color='r', label='Valeur exacte de l\'intégrale')
plt.xlabel('N : Nombre de subdivisions')
plt.ylabel('Valeur de l\'intégrale par les différentes méthodes d\'intégration')
plt.title('Comparaison d\'évolution de la valeur de l\'intégrale par les différentes méthodes d\'intégration')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.show()

