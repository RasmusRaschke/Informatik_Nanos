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
    """cheap method to find zeros in wavefunction by sign change

        Parameters
        ----------
        u_0 : float
            first boundary value
        u_1 : float
            second fixed point for numerov algorithm
        v : float list
            potential well defined as 1-D array
        counter : int
            uses counter (=len(v)) to make sure that index range is conserved when using a variable amount of wells
        x_min : float
            start of x-axis
        x_max : float
            end of x-axis
        N : int
            number of steps
        e_min : float
            starting point for energy search
        e_max : float
            ending point for energy search
        accuracy: float
            accuracy of search for zeros

        Returns
        -------
        eigen : float list
            rough eigenenergies of the wave function
        """
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
    return eigen


def newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, eigen_estimate, accuracy, differential):
    """algorithm to increase the accuracy of found zeros

        Parameters
        ----------
        u_0 : float
            first boundary value
        u_1 : float
            second fixed point for numerov algorithm
        counter : int
            uses counter (=len(v)) to make sure that index range is conserved when using a variable amount of wells
        x_min : float
            start of x-axis
        x_max : float
            end of x-axis
        N : int
            number of steps
        v : float list
            potential well defined as 1-D array
        max_bound : int
            maximal amount of iterations if control doesn't get low enough
        eigen_estimate : float
            eigenenergy from the inaccurate method
        accuracy : float
            convergence criterium for zeros
        differential : float
            first step to calculate the new energy

        Returns
        -------
        e_new : float list
            nearly exact eigenenergies of the wave function
        """
    e = eigen_estimate
    de = differential
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

def correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy, differential, eigen_estimate):
    """algorithm to increase the accuracy of found zeros

        Parameters
        ----------
        u_0 : float
            first boundary value
        u_1 : float
            second fixed point for numerov algorithm
        counter : int
            uses counter (=len(v)) to make sure that index range is conserved when using a variable amount of wells
        x_min : float
            start of x-axis
        x_max : float
            end of x-axis
        N : int
            number of steps
        v : float list
            potential well defined as 1-D array
        max_bound : int
            maximal amount of iterations if control doesn't get low enough
        accuracy : float
            convergence criterium for zeros
        differential : float
            first step to calculate the new energy
        eigen_estimate : float
            eigenenergy from the inaccurate method

        Returns
        -------
        eigen_corrected : float list
            nearly exact eigenenergies of the wave function
        """
    eigen_corrected = []
    for i in range(len(eigen_estimate)):
        e = newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, eigen_estimate[i], accuracy, differential)
        eigen_corrected.append(e)
    return eigen_corrected


def truncate(u, factor):
    maximum = max(u)
    limit = maximum / factor
    for i in range(len(u)):
        if abs(u[i]) < limit:
            u[i] = 0
    return u


def calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max, accuracy, max_bound, start):
    x, v, counter = periodic(v_0, x_0, L, N, x_min, x_max, wells)
    eigen_simple = find_zeros(u_0, u_1, v, counter, x_min, x_max, N, e_min, e_max, accuracy)
    eigen_correct = correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy, start, eigen_simple)
    return eigen_correct, x, v, counter


def plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, N):
    fig, ax = plt.subplots()
    ax.plot(dpi=300)
    ax.plot(x, v)
    ax2 = ax.twinx()
    for i in range(len(eigen_correct)):
        K = get_k(v, eigen_correct[i])
        wave = numerov(u_0, u_1, K, counter, x_min, x_max, N)
        ax2.plot(x, norm(wave), label=r"$\psi_{%i}$"%(i+1))
    plt.grid(True)
    ax.set_xlabel(r'x in $[nm]$')
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red')
    plt.legend()
    plt.show()


def plot_bands():
    pass
########################################################################################################################
#Variablen hier
v_0 = 0.3
x_0 = 2
L = 3
N = 1000
x_min = 0
x_max = 10
wells = 50
E = 1.
u_0 = .0
u_1 = .01
max_bound = 1000000
accuracy = 0.00001
start = 0.01
E_0 = -1
E_max = 0
e_min = -2
e_max = 0
########################################################################################################################
#Berechnung
eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,
                                                     accuracy, max_bound, start)
plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, N)




