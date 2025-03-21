"""

@author: Alexandre Nuyt

"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jn
from scipy.optimize import newton

from sage.symbolic.operators import *
from sage.functions.special import *

xi, k, h, n, a = var('xi k h n a')
J_k = bessel_J(abs(k), xi)

equation = xi*(-J_k + k/xi * bessel_J(abs(k+1), xi)) +
 k*(1 + (xi^2 * h^2)/(n^2 * pi^2 * a^2))^(1/2) * bessel_J(abs(k), xi) == 0

# Dérivée par rapport à xi
derivative_equation = diff(equation, xi)

show(derivative_equation)

# Newton-Raphson   

# 1. Dérivée Numérique

# Paramètres
k = 1
h = 1
n = 1
a = 1

xi_values = np.linspace(0.001, 25, 10000)

def equation(xi, k, h, n, a):
    eq = xi * (-jn(np.abs(k+1), xi) + k/xi * jn(np.abs(k), xi)) +
    k * np.sqrt((1 + (xi**2 * h**2) / (n**2 * np.pi**2 * a**2))) * 
    jn(np.abs(k), xi)
    return eq
    
def d_equation(xi, k, h, n, a):
    h = 1e-8
    return (equation(xi + h, k, h, n, a) - equation(xi, k, h, n, a)) / h
    
"""
# 2. Dérivée analytique

def derivee_equation(xi, k, h, n, a):
    d_eq = -1/2*k*np.sqrt(h^2*xi^2/(np.pi^2*a^2*n^2) + 1)*(jn(np.abs(k+1), xi) 
    - jn(np.abs(k-1), xi)) - 1/2*xi*(k*(jn(np.abs(k + 1) + 1, xi) - 
    jn(np.abs(k + 1) - 1, xi))/xi + 2*k*jn(np.abs(k + 1), xi)/xi^2 - 
    jn(np.abs(k) + 1, xi) + jn(np.abs(k) - 1, xi)) +
    k*jn(np.abs(k + 1), xi)/xi + h^2*k*xi*jn(np.abs(k), xi)/
    (np.pi^2*a^2*n^2*np.sqrt(h^2*xi^2/(np.pi^2*a^2*n^2) + 1)) - 
    jn(np.abs(k), xi)
    return d_eq  
"""

def newton_raphson(k, h, n, a, nb_racines=7):
    racines = []
    xi_ini = 0.1
    
    for _ in range(nb_racines):
        xi_0_NR = newton(equation, xi_ini, args=(k, h, n, a),
                             fprime=d_equation, tol=1e-8, maxiter=1000)
        racines.append(xi_0_NR)
        xi_ini = xi_0_NR + 3.0  
        # Utiliser la solution trouvée comme nouvelle supposition initiale
    
    return racines

# Fonction
equation_values = equation(xi_values, k, h, n, a)

# Racines
racines = newton_raphson(k, h, n, a, nb_racines=7)

# Tracer la fonction et les racines
plt.figure(figsize=(6, 4))
plt.plot(xi_values, equation_values, label='Equation transcendante')
plt.scatter(racines, [0]*len(racines), color='red', label='Racines trouvées',
            marker='*')
plt.title(f'Racines de l\'équation pour $k={k}$, $h={h}$, $n={n}$, $a={a:.1f}$',
          fontsize=10)
plt.xlabel('$\\xi$')
plt.ylabel('Fonction')
for xi in racines:
    plt.text(xi, 0, f'$\\xi={xi:.3f}$', fontsize=8, ha='right', va='bottom')
plt.legend()
plt.grid(True)
plt.show()

# Calcul des valeurs porpres
lambda_values = []
for xi_0 in racines : 
    lambda_nmk = 2 / np.sqrt(1 + ((xi_0)**2 * h**2) / (n**2 * np.pi**2 * a**2))
    lambda_values.append(lambda_nmk)
    
# Tracer le graphe des valeurs de lambda en fonction des racines
plt.figure(figsize=(6, 4))
plt.scatter(racines, lambda_values, color='blue', label='$\\lambda_{nmk}$')
plt.title(f'Valeurs de $\\lambda$ pour $k={k}$, $h={h}$, $n={n}$, $a={a:.1f}$',
          fontsize=10)
plt.xlabel('$\\xi$')
plt.ylabel('$\\lambda_{nmk}$')
for xi, lambd in zip(racines, lambda_values):
    plt.text(xi, lambd, f'(${xi:.2f}$, ${lambd:.2f}$)', fontsize=8, ha='right',
             va='bottom')
plt.legend()
plt.grid(True)
plt.show()
