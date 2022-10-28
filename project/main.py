import numpy as np
import matplotlib.pyplot as plt
import math as m
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
        k.append((v[i] - E) * 26.27)
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
    h = (x_max - x_min) / N
    for i in range(2, counter):
        u.append((2 * u[i - 1] * (1 - (5 / 12) * h**2 * K[i - 1]) - u[i - 2] * (1 + h**2 / 12 * K[i - 2])) / (1 + 1 / 12 * h**2 * K[i]))
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
        sum_ += u[i]**2
        if sum != 0:
            u_norm.append(u[i] / np.sqrt(sum_))
        else:
            u_norm.append(0)
    return u_norm


def newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, energy, accuracy, start):
    e = energy
    de = start
    for i in range(max_bound):
        K = get_k(v, e)
        psi_1 = numerov(u_0, u_1, K, counter, x_min, x_max, N)
        psi_max = psi_1[N-1]
        K = get_k(v, e + de)
        psi_2 = numerov(u_0, u_1, K, counter, x_min, x_max, N)
        psi_max_de = psi_2[N-1]
        de_new = - psi_max / ((psi_max_de - psi_max) / de)
        e_new = e + de_new
        de = de_new
        control = e_new - e
        e = e_new
        if control < accuracy:
            break
    return e_new


########################################################################################################################
#Variablen hier
v_0 = 5
x_0 = 1
L = 2
N = 1000
x_min = 0
x_max = 10
wells = 5
E = 0.1
u_0 = .0
u_1 = .01
max_bound = 1000000
accuracy = 0.01
start = 0.1
########################################################################################################################
#Berechnung
x, v, counter = periodic(v_0, x_0, L, N, x_min, x_max, wells)
eigen = newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, E, accuracy, start)
K = get_k(v, eigen)
psi = numerov(u_0, u_1, K, counter, x_min, x_max, N)
print(psi)
psi_norm = norm(psi)
print(eigen)
########################################################################################################################
#Plot
fig, ax = plt.subplots()
ax.plot(x, v)
ax.plot(x, get_k(v, -E))
ax.plot(x, psi_norm)
plt.show()

