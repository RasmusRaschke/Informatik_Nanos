import numpy as np
import matplotlib.pyplot as plt
########################################################################################################################
def x_axis(N, x_min, x_max):
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

def get_k(N, x_min, x_max, v, E):
    k = []
    h = (x_max - x_min) / N
    for i in range(N):
        k.append((v[i] - E) * 26.27)
    return k

def numerov(u_0, u_1, K, N, x_min, x_max):
    u = [u_0, u_1]
    h = (x_max - x_min) / N
    for i in range(2, N):
        u.append((2 * u[i - 1] * (1 - (5 / 12) * h**2 * K) - u[i - 2] * (1 + h**2 / 12 * K)) / (1 + 1 / 12 * h**2 * K))
    return u

def norm()
########################################################################################################################


