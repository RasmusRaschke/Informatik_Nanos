import numpy as np
import matplotlib.pyplot as plt
import math as m
#print(plt.style.available)
#plt.style.use('ggplot')
'''
to do: change N_inside and N_outside to N_ges (easier for Numerov), complete Numerov in C++
usage: tbw
'''

def euler(u_0, z_0, K, N, x_min, x_max):
    u = []
    u.append(u_0)
    z = []
    z.append(z_0)
    h = (x_max - x_min) / N
    for i in range(1, N):
        u.append(z[i - 1] * h + u[i - 1])
        z.append(z[i-1] - h * K * u[i-1])
    return u


def numerov(u_0, u_1, K, N, x_min, x_max):
    u = [u_0, u_1]
    h = (x_max - x_min) / N
    for i in range(2, N):
        u.append((2 * u[i - 1] * (1 - (5 / 12) * h**2 * K) - u[i - 2] * (1 + h**2 / 12 * K)) / (1 + 1 / 12 * h**2 * K))
    return u

def x_axis(N_outside, N_inside, x_min, x_max):
    N = 2 * N_outside + N_inside
    h = (x_max - x_min) / N
    x = []
    for i in range(N):
        x.append(i * h)
    return x

def potential(v_0, x_0, L, N_outside, N_inside, x_min, x_max):
    width = x_max - x_min
    left_width = L / 2
    h_left = left_width / N_outside
    h_inside = L / N_inside
    v = np.zeros(N_outside)
    for i in range(N_inside):
        v = np.append(v, - v_0)
    v = np.append(v, np.zeros(N_outside))
    return v

def get_k(N_outside, N_inside, x_min, x_max, v, E):
    k = []
    N = 2 * N_outside + N_inside
    h = (x_max - x_min) / N
    for i in range(N):
        k.append((v[i] - E) * 25.27)
    return k

def wavefunction


v, x = potential(1, 2, 5, 100, 200, 0, 10), x_axis(100, 200, 0, 10)
k = get_k(100, 200, 0, 10, v, - 0.0034)
fig, ax = plt.subplots()
ax.plot(x, v, label="Potential")
ax.plot(x, k, label="K-Liste")
ax.set_xlabel(r"$l$ in $[nm]$")
ax.set_ylabel(r"$V$ in $[eV]$")
ax.legend()
ax.grid(True)
plt.show()



