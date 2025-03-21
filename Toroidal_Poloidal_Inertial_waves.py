"""

@author: Alexandre Nuyt

"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import time

class Inertial_Waves_Calculator:
    def __init__(self, n=1000, s_min=0, s_max=2, Omega=1):
        self.n = n
        self.s_min = s_min
        self.s_max = s_max
        self.Omega = Omega
        self.h = (s_max - s_min) / n
        self.s = np.linspace(s_min, s_max, 2*n)
        self.k = None
        self.m = None
        self.A = None
        self.B = None
        self.eigenvalues_NS_sorted = None
        self.eigenvectors_NS = None
        self.omega = None
        self.index_omega = None
        
    def input_parameters(self):
        self.k = int(input("Insérer une valeur de k : "))
        self.m = int(input("Insérer une valeur de m : "))
        
    def matrix_A1(self):
        return np.zeros((self.n, self.n))
    def matrix_A2(self):
        A2 = np.zeros((self.n, self.n), dtype=complex)
        for i in range(self.n):
            if i == self.n - 1:
                A2[i, i] = 1
            else:
                A2[i, i] = (-2 * self.Omega * self.k * (self.m ** 2) * 
                            self.s[i] + (self.s[i] ** 2) * (self.k ** 3)) * 1j
        return A2
    def matrix_A3(self):
        A3 = np.zeros((self.n, self.n))
        for i in range(self.n):
            A3[i, i] = self.s[i] ** 2
        return (2 * self.Omega * (self.k ** 2)) * A3
    def matrix_A4(self):
        A4 = np.zeros((self.n, self.n))
        for i in range(self.n):
            if i == 0:
                A4[i, i] = -(1 + (self.s[i]/(2 * self.h)))
                A4[i, i + 1] = self.s[i] / (2 * self.h)
            elif i == self.n - 1:
                A4[i, i] = 0
            else:
                A4[i, i - 1] = -self.s[i] / (2 * self.h)
                A4[i, i] = -1
                A4[i, i + 1] = self.s[i] / (2 * self.h)
        return (2 * self.Omega * self.k * self.m) * A4

    def matrix_A(self):
        A1 = self.matrix_A1()
        A2 = self.matrix_A2()
        A3 = self.matrix_A3()
        A4 = self.matrix_A4()

        A1_2 = np.block([[A1, A2]])
        A3_4 = np.block([[A3, A4]])
        self.A = np.block([[A1_2], [A3_4]])

    def matrix_B1(self):
        B1 = np.zeros((self.n, self.n), dtype=complex)
        for i in range(self.n):
            B1[i, i] = ((self.m ** 2) + (self.k ** 2) * (self.s[i] ** 2)) * 1j
        return B1
    def matrix_B2(self):
        B2 = np.zeros((self.n, self.n), dtype=complex)
        for i in range(self.n):
            if i == self.n - 1:
                B2[i, i] = 1
            else:
                B2[i, i] = -2j
        return B2
    def matrix_B3(self):
        B3 = np.zeros((self.n, self.n))
        for i in range(self.n):
            if i == 0:
                B3[i, i] = -(1 + (self.s[i] / (2 * self.h)))
                B3[i, i + 1] = self.s[i] / (2 * self.h)
            elif i == self.n - 1:
                B3[i, i] = -(1 + (self.s[i] /( 2 * self.h)))
                B3[i, i - 1] = self.s[i] / (2 * self.h)
            else:
                B3[i, i - 1] = -self.s[i] / (2 * self.h)
                B3[i, i] = -1
                B3[i, i + 1] = self.s[i] / (2 * self.h)
        return self.m * B3
    def matrix_B4_1(self):
        B4_1 = np.zeros((self.n, self.n))
        for i in range(self.n):
            if i == self.n - 1:
                B4_1[i, i] = 1
            else:
                B4_1[i, i] = self.s[i] * self.k * (self.m ** 2) + 
                (self.s[i] ** 2) + (self.k ** 3) + self.k
        return B4_1
    def matrix_B4_2(self):
        B4_2 = np.zeros((self.n, self.n))
        for i in range(self.n):
            if i == 0:
                B4_2[i, i] = 2 * (self.s[i] ** 2) + self.h * self.s[i]
                B4_2[i, i + 1] = -4 * (self.s[i] ** 2) + self.h * self.s[i]
                B4_2[i, i + 2] = 2 * (self.s[i] ** 2)
            elif i == self.n - 1:
                B4_2[i, i] = 1
            else:
                B4_2[i, i - 1] = 2 * (self.s[i] ** 2) + self.h * self.s[i]
                B4_2[i, i] = -4 * (self.s[i] ** 2)
                B4_2[i, i + 1] = 2 * (self.s[i] ** 2) + self.h * self.s[i]
        return (self.k / (2 * (self.h ** 2))) * B4_2

    def matrix_B(self):
        B1 = self.matrix_B1()
        B2 = self.matrix_B2()
        B3 = self.matrix_B3()
        B4 = self.matrix_B4_1() + self.matrix_B4_2()

        B1_2 = np.block([[B1, B2]])
        B3_4 = np.block([[B3, B4]])
        self.B = np.block([[B1_2], [B3_4]])

    def calculate_eigenvalues_and_vectors(self):
        self.eigenvalues_NS, self.eigenvectors_NS = linalg.eig(self.A, self.B)
        self.eigenvalues_NS_sorted = sorted(self.eigenvalues_NS)
        index_eigenvalues_NS_sorted = np.argsort(self.eigenvalues_NS)
        self.index_omega = index_eigenvalues_NS_sorted[0]
        self.omega = np.real(np.abs(self.eigenvalues_NS_sorted[0]))
        print('Pour n =', self.n)
        print('omega_0', self.omega)
        print('omega_1', np.real(np.abs(self.eigenvalues_NS_sorted[1])))
        print('omega_2', np.real(np.abs(self.eigenvalues_NS_sorted[2])))
        print('omega_3', np.real(np.abs(self.eigenvalues_NS_sorted[3])))

        plt.figure(figsize=(10, 8))
        plt.rcParams.update({'font.family': 'arial'})
        plt.subplot(3, 2, 3)
        vep_Toro = (self.s[:len(self.s)//2]**2)*
        self.eigenvectors_NS[:len(self.s) // 2, self.index_omega]
        plt.plot(self.s[:len(self.s)//2], vep_Toro, label="Vep Toroïdal",
                 color='blue')
        plt.xlabel('Position $s$')
        plt.ylabel('Amplitude')
        plt.title('Composante Toroïdale')
        plt.grid()
        plt.legend()
        plt.subplot(3, 2, 4)
        vep_Polo = (self.s[:len(self.s)//2]**2)*
        self.eigenvectors_NS[len(self.s) // 2:, self.index_omega]
        plt.plot(self.s[:len(self.s)//2], vep_Polo, label="Vep Poloïdal",
                 color='green')
        plt.xlabel('Position $s$')
        plt.ylabel('Amplitude')
        plt.title('Composante Toroïdale')
        plt.grid()
        plt.legend()
        plt.tight_layout()
        plt.show()

    def calculate_velocity_components(self):
        z_min = 0
        z_max = 10
        t_0 = 0
        t_max = (2*np.pi) / self.omega

        phi = np.linspace(0, 2 * np.pi, self.n)
        z = np.linspace(z_min, z_max, self.n)
        t = np.linspace(t_0, t_max, self.n)

        v_s = []
        v_phi = []
        v_z = []

        vep_Toro = (self.s[:len(self.s)//2]**2)*
        self.eigenvectors_NS[:len(self.s) // 2, self.index_omega]
        vep_Polo = (self.s[:len(self.s)//2]**2)*
        self.eigenvectors_NS[len(self.s) // 2:, self.index_omega]

        for i in range(self.n):
            if self.s[i] != 0:
                exp_term = np.exp(1j * (self.k * z[i] + self.m * phi[i] +
                                        self.omega * t[i]))
                if i == 0:
                    v_s.append((self.s[i] ** 2) * exp_term * 
                    ((-1 / (self.s[i] ** 2)) * (-self.m ** 2) * vep_Polo[i] - 
                    (-self.k ** 2) * vep_Polo[i]))
                    v_phi.append((self.s[i] ** 2) * exp_term *
                    (1j * self.k * vep_Toro[i] + (1j * self.m / self.s[i]) *
                    (1 / self.h) * (vep_Polo[i + 1] - vep_Polo[i]) - 
                    (1j * self.m / self.s[i] ** 2) * vep_Polo[i]))
                    v_z.append(self.s[i] ** 2 * exp_term * 
                    ((1j / self.s[i]) * (self.k * vep_Polo[i] - 
                    self.m * vep_Toro[i]) + (1j * self.k / self.h) * 
                    (vep_Polo[i + 1] - vep_Polo[i])))
                elif i == self.n - 1:
                    v_s.append(self.s[i] ** 2 * exp_term * 
                    ((-1 / (self.s[i] ** 2)) * (-self.m ** 2) * 
                    vep_Polo[i] - (-self.k ** 2) * vep_Polo[i]))
                    v_phi.append(self.s[i] ** 2 * exp_term * 
                    (1j * self.k * vep_Toro[i] + (1j * self.m / self.s[i]) * 
                    (1 / self.h) * (vep_Polo[i] - vep_Polo[i - 1]) - 
                    (1j * self.m / self.s[i] ** 2) * vep_Polo[i]))
                    v_z.append(self.s[i] ** 2 * exp_term * 
                    ((1j / self.s[i]) * (self.k * vep_Polo[i] - 
                                         self.m * vep_Toro[i]) + 
                    (1j * self.k / self.h) * (vep_Polo[i] - vep_Polo[i - 1])))
                else:
                    v_s.append(self.s[i] ** 2 * exp_term * 
                    ((-1 / (self.s[i] ** 2)) * (-self.m ** 2) * vep_Polo[i] - 
                    (-self.k ** 2) * vep_Polo[i]))
                    v_phi.append(self.s[i] ** 2 * exp_term * 
                    (1j * self.k * vep_Toro[i] + (1j * self.m / self.s[i]) * 
                    (1 / (2 * self.h)) * (vep_Polo[i + 1] - vep_Polo[i - 1]) - 
                    (1j * self.m / self.s[i] ** 2) * vep_Polo[i]))
                    v_z.append(self.s[i] ** 2 * exp_term * ((1j / self.s[i]) * 
                    (self.k * vep_Polo[i] - self.m * vep_Toro[i]) + 
                    (1j * self.k / (2 * self.h)) * (vep_Polo[i + 1] -
                    vep_Polo[i - 1])))

        return v_s, v_phi, v_z

    def plot_velocity_components(self, v_s, v_phi, v_z):
        plt.rcParams.update({'font.family': 'arial'})
        fig, axs = plt.subplots(3, 1, sharex=True)
        fig.subplots_adjust(hspace=0.5)

        axs[0].plot(self.s[:len(self.s) // 2 - 2], v_s[:-1])
        axs[0].set_ylim(np.nanmin(v_s), np.nanmax(v_s))
        axs[0].set_title('$\sim v_s$ en fonction de s')

        axs[1].plot(self.s[:len(self.s) // 2 - 2], v_phi[:-1])
        axs[1].set_ylim(np.nanmin(v_phi), np.nanmax(v_phi))
        axs[1].set_title('$\sim v_\phi$ en fonction de s')

        axs[2].plot(self.s[:len(self.s) // 2 - 2], v_z[:-1])
        axs[2].set_ylim(np.nanmin(v_z), np.nanmax(v_z))
        axs[2].set_title('$\sim v_z$ en fonction de s')
        axs[2].set_xlabel('Position $s$')

        plt.show()

    def plot_v_phi_propagation(self, v_phi):
        z = np.linspace(0, 10, self.n)
        plt.rcParams.update({'font.family': 'arial'})
        plt.figure(figsize=(8, 6))
        plt.scatter(self.s[:len(self.s) // 2 - 1], z[:len(z) - 1], c=v_phi[:],
                    cmap='viridis')
        plt.colorbar(label='Amplitude des oscillations de $v_{\phi}$')
        plt.xlabel('Position $s$')
        plt.ylabel('$z$')
        plt.title('Amplitude des oscillations de $v_{\phi}$')
        plt.grid(False)
        plt.show()

    def run(self):
        start_time = time.time()

        self.input_parameters()
        self.matrix_A()
        self.matrix_B()
        self.calculate_eigenvalues_and_vectors()

        v_s, v_phi, v_z = self.calculate_velocity_components()
        self.plot_velocity_components(v_s, v_phi, v_z)
        self.plot_v_phi_propagation(v_phi)

        end_time = time.time()
        execution_time = end_time - start_time
        print("Temps d'exécution:", execution_time, "secondes")

if __name__ == "__main__":
    calculator = Inertial_Waves_Calculator()
    calculator.run()

