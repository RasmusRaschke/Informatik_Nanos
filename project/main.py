import numpy as np
import matplotlib.pyplot as plt
import math as m
import copy
########################################################################################################################
def periodic(v_0, x_0, L, N, x_min, x_max, wells=1):
    """Get equally spaced potential wells of defined width and depth; calculate correspondent x-axis in nm

    Parameters
    ----------
    v_0 : float
        depth of potential wells in eV; inversed wells (walls) aren't supported right now
    x_0 : float
        point where the first well begins
    L : float
        width of the well
    N : int
        steps for resolution
    x_min : float
        start of x-axis
    x_max : float
        end of x-axis
    wells : int, optional
        number of equally spaced wells, default = 1

    Returns
    -------
    v : float list
        1-D potential array
    x : float list
        x-axis in nm
    counter : int
        returns len of v and x
    """
    #error handling: set walls of well to 0
    if v_0 < 0:
        v_0 = - v_0
    v = []
    x = []
    h = (x_max - x_min) / N
    spacing_steps = m.ceil((x_0 - x_min) / h)
    well_steps = m.ceil(L / h)
    counter = 0
    for i in range(wells):
        for j in range(spacing_steps):
            v.append(0)
            counter += 1
        for k in range(well_steps):
            v.append(- v_0)
            counter += 1
    for i in range(spacing_steps):
        v.append(0)
        counter += 1
    for i in range(counter):
        x.append(i * h)
    return x, v, counter

def get_k(v, E):
    """Get k-list

    Parameters
    ----------
    v : float list
        potential well defined as 1-D array
    E : float
        energy
    Returns
    -------
    k : float list
        1-D energy array
    """
    k = []
    for i in range(len(v)):
        k.append(- (v[i] - E) * 26.27)
    return k

def numerov(u_0, u_1, K, counter, x_min, x_max, N):
    """Apply numerov algorithm to solve separated SEQ

    Parameters
    ----------
    u_0 : float
        first boundary value
    u_1 : float
        second fixed point for numerov algorithm
    K : float
        constant in SQE
    counter : int
        uses counter (=len(v)) to make sure that index range is conserved when using a variable amount of wells
    x_min : float
        start of x-axis
    x_max : float
        end of x-axis
    N : int
        number of steps

    Returns
    -------
    u : float
        wave function as 1-D array calculated with numerov

    """
    u = [u_0, u_1]
    h = (x_max - x_min) / counter
    for i in range(2, counter):
        u.append(((2 * u[i - 1] * (1 - (5 / 12) * h**2 * K[i - 1])) - u[i - 2] * (1 + h**2 * K[i - 2] / 12)) / (1 + 1 / 12 * h**2 * K[i]))
    return u

def norm(u):
    """returns normed wave function

        Parameters
        ----------
        u : list
            wave function from numerov algorithm

        Returns
        -------
        u_norm : list
            normed wave function

        """
    sum_ = 0
    u_norm = []
    for i in range(len(u)):
        sum_ = u[i]**2 + sum_
    for i in range(len(u)):
        u_norm.append(u[i] / np.sqrt(sum_))
    return u_norm

def find_zeros(u_0, u_1, v, counter, x_min, x_max, N, e_min, e_max, accuracy):
    energies = []
    eigen = []
    last = []
    n = abs(round((e_max - e_min) / accuracy))
    for i in range(n):
        energy = e_min + i * accuracy
        K = get_k(v, energy)
        psi = numerov(u_0, u_1, K, counter, x_min, x_max, N)
        last.append(psi[-1])
        energies.append(energy)
    for i in range(1, len(energies)):
        if np.sign(last[i-1]) != np.sign(last[i]):
            energy = (i * accuracy - 0.5 * accuracy) + e_min
            eigen.append(energy)
    return(eigen)


def newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, energy, accuracy, start):
    e = energy
    de = start
    for i in range(max_bound):
        K = get_k(v, e)
        psi_1 = numerov(u_0, u_1, K, counter, x_min, x_max, N)
        psi_max = psi_1[-1]
        K = get_k(v, e + de)
        psi_2 = numerov(u_0, u_1, K, counter, x_min, x_max, N)
        psi_max_de = psi_2[-1]
        de_new = - psi_max / ((psi_max_de - psi_max) / de)
        e_new = e + de_new
        de = de_new
        control = e_new - e
        e = copy.copy(e_new)
        if control < accuracy:
            break
    print(e_new)
    return e_new

def correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy, start, zeros):
    eigen_correct = []
    for i in range(len(zeros)):
        e = newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, zeros[i], accuracy, start)
        eigen_correct.append(e)
    return eigen_correct

########################################################################################################################
#Variablen hier
v_0 = 1
x_0 = 2
L = 5
N = 1000
x_min = 0
x_max = 10
wells = 2
E = 1.
u_0 = .0
u_1 = .01
max_bound = 1000000
accuracy = 0.001
start = 0.01
E_0 = -1
E_max = 0
e_min = -2
e_max = 0
########################################################################################################################
#Berechnung
x, v, counter = periodic(v_0, x_0, L, N, x_min, x_max, wells)
print(counter)
eigen_simple = find_zeros(u_0, u_1, v, counter, x_min, x_max, N, e_min, e_max, accuracy)
eigen_correct = correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy, start, eigen_simple)
print("Einfacher Wert:", eigen_simple)
print("Genauer Wert:", eigen_correct)
K1 = get_k(v, eigen_correct[0])
psi1 = numerov(u_0, u_1, K1, counter, x_min, x_max, N)
K2 = get_k(v, eigen_correct[1])
psi2 = numerov(u_0, u_1, K2, counter, x_min, x_max, N)
K3 = get_k(v, eigen_correct[2])
psi3 = numerov(u_0, u_1, K3, counter, x_min, x_max, N)
plt.plot(x, norm(psi1))
plt.plot(x, norm(psi2))
plt.plot(x, norm(psi3))
plt.show()
########################################################################################################################
#Plot
fig, ax = plt.subplots()
ax.plot(dpi=300)
ax.plot(x, v)
ax.plot(x, norm(psi1))
plt.show()

