"""

@author: Alexandre Nuyt

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jn_zeros

A = 1.0
L = 1.0
m_values = [0, 1, 2, 3, 4,5]
s_values = np.linspace(0, L, 1000)
def J(s, m):
    return A * jv(m, s)
"""
Les valeurs propres sont données par les zéros des fonctions de Bessel
""" 
for m in m_values:
    # Racines Bessel
    zeros = jn_zeros(m, L)
    print("Vap pour m ="f'{m}')
    for zero in zeros:
        print(round(zero, 5)) #valeurs arrondies à cinq décimales
    
    plt.plot(s_values, J(s_values*zeros[0], m), label=f'$J_{m!r}$')
plt.xlabel('s')
plt.ylabel('$J_{m}(ks)$')
plt.title('Fonctions de Bessel')
plt.grid(True)
plt.legend()
plt.show()

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

A = 1.0
n = 500
s_min = 0
s_max = 1.0
s = np.linspace(s_min, s_max, n)
h = s_max/n
m_values = [0,1,2,3,4,5]

""" 
Ici on impose comme condition initiale que df(s)/ds en s=0 vaut 0 et que 
f(s_n) = 0 en s = 1
"""
def matrix_A_1(n, h, m):
    if m % 2 == 0:
        A_1 = np.zeros((n,n))            
        for i in range(n):
            if i == 0:
                A_1[i, i] = -2
                A_1[i,i+1] = 2
            elif i == n - 1:
                A_1[i, i] = 2*h
            else:
                A_1[i, i-1] = -1*s[i]
                A_1[i, i] = 0
                A_1[i, i+1] = 1*s[i]
    else: 
        A_1 = np.zeros((n,n))
        for i in range(n):
            if i == 0:
                A_1[i, i] = 2*h
            elif i == n - 1:
                A_1[i, i] = 2*h
            else:
                A_1[i, i-1] = -1*s[i]
                A_1[i, i] = 0
                A_1[i, i+1] = 1*s[i]
    return (1 / (2*h)) * A_1

"""
On impose les premières et dernières lignes à zéro afin de faire valoir les
 conditions aux bords
"""
def matrix_A_2(n, h):
    A_2 = np.zeros((n,n))
    for i in range(n):
        if i == 0:
            A_2[i, i] = 0
        elif i == n - 1:
            A_2[i, i] = 0
        else:
            A_2[i, i-1] = 1*s[i]**2
            A_2[i, i] = -2*s[i]**2
            A_2[i, i+1] = 1*s[i]**2
    return (1 / h**2) * A_2

def matrix_A_3(n,h,m):
    A_3 = np.zeros((n,n))
    for i in range(n):
        if i == 0:
            A_3[i, i] = 0
        elif i == n - 1:
            A_3[i, i] = 0
        else:
            A_3[i,i] = 1
    return (m**2)*A_3

"""
On impose ici les valeurs que doit prendre la fonction en s=0,1
"""
import numpy as np

def matrix_B(n, h, m):
    if m % 2 == 0:
        B = np.zeros((n, n))
        for i in range(n):
            if i == 0:
                B[i, i] = 0
                B[i,i+1] = 0
            elif i == 1:
                B[i, i] = 0
            elif i == n-1:
                B[i, i] = 0
            else:
                B[i, i] = -s[i]**2
    else:
        B = np.zeros((n, n))
        for i in range(n):
            if i == 0:
                B[i, i] = 0
            elif i == 1:
                B[i, i] = 0
            elif i == n-1:
                B[i, i] = 0
            else:
                B[i, i] = -s[i]**2
    return B
"""
for m in m_values:
    A = matrix_A_1(n,h,m) + matrix_A_2(n,h) - matrix_A_3(n,h,m)
    B = matrix_B(n,h,m)

#print(A)
#print(B)

    eigenvalues, eigenvectors = linalg.eig(A, B)
    eigenvalues_sorted = sorted(eigenvalues)

    index_eigenvalues_sorted = np.argsort(eigenvalues) # renvoie les indices 
    des vap triées par ordre croissant

    plt.figure(figsize=(6, 4))
    #for i in range(0,1): 
    plt.plot(s, np.abs(eigenvectors[:, index_eigenvalues_sorted[0]]) , 
             label=f"Vecteur propre pour m = {m}")
    print(f"Pour m = {m}, 
          la première valeur propre est : {np.sqrt((eigenvalues_sorted[i]))})")

    plt.xlabel('Position')
    plt.ylabel('Amplitude')
    plt.title('Vecteurs propres numériques triés')
    plt.legend(loc='best')
    plt.grid(True)

    plt.show()
"""
first_eigenvalues = []
first_eigenvectors = []

for m in m_values:
    A = matrix_A_1(n, h, m) + matrix_A_2(n, h) - matrix_A_3(n, h, m)
    B = matrix_B(n, h, m)

    eigenvalues, eigenvectors = linalg.eig(A, B)
    eigenvalues_sorted = np.sort(eigenvalues)
    index_eigenvalues_sorted = np.argsort(eigenvalues)

    first_eigenvalues.append(round(np.real(np.sqrt(eigenvalues_sorted[0])),5))
    first_eigenvectors.append(eigenvectors[:, index_eigenvalues_sorted[0]])

#plt.figure(figsize=(8, 6))

for i, m in enumerate(m_values):
    print(f"Pour m = {m}, la première valeur propre est : {first_eigenvalues[i]}")
    
"""
plt.plot(s, np.abs(first_eigenvectors[i]), label=f"m = {m}, 
         $\omega$ = {first_eigenvalues[i]}")

plt.xlabel('Position')
plt.ylabel('Amplitude')
plt.title('Premières valeurs propres pour différentes valeurs de m')
plt.legend(loc='best')
plt.grid(True)

plt.show()
"""
