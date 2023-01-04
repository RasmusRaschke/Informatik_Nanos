import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import math as m
import copy
import tkinter as tk
from tkinter import ttk
########################################################################################################################

# create a window that will be needed to show the results and change the variables

window = tk.Tk()


on = 1

# set variables in case the User don´t use the advance Options

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

def numerov(u_0, u_1, K, counter, x_min, x_max):
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
        psi = numerov(u_0, u_1, K, counter, x_min, x_max)
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
    c = 0
    for i in range(max_bound):
        K = get_k(v, e)
        psi_1 = numerov(u_0, u_1, K, counter, x_min, x_max)
        psi_max = psi_1[-1]
        K = get_k(v, e + de)
        psi_2 = numerov(u_0, u_1, K, counter, x_min, x_max)
        psi_max_de = psi_2[-1]
        de_new = - psi_max / ((psi_max_de - psi_max) / de)
        e_new = e + de_new
        de = de_new
        control = e_new - e
        e = copy.copy(e_new)
        c += 1
        if control < accuracy:
            break
    print(e_new)
    print("Schritte:", c)
    return e_new

def correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy, start, eigen_estimate):
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
        start : float
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
        e = newton_raphson(u_0, u_1, counter, x_min, x_max, N, v, max_bound, eigen_estimate[i], accuracy, start)
        eigen_corrected.append(e)
    return eigen_corrected


def calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max, accuracy, max_bound, start):
    """Calculate accurate eigenvalues with NR-algorithm

        Parameters
        ----------
        v_0 : float
            depth of potential wells in eV; inverse wells (walls) aren't supported right now
        x_0 : float
            point where the first well begins
        u_0 : float
            first boundary value
        u_1 : float
            second fixed point for numerov algorithm
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
        e_min : float
            starting point for energy search
        e_max : float
            ending point for energy search
        accuracy : float
            accuracy of search for zeros
        max_bound : int
            maximal amount of iterations if control doesn't get low enough
        start : float
            first step to calculate the new energy

        Returns
        -------
        eigen_correct : float list
            list of eigenenergies in [eV]
        x : float list
            x-axis in nm
        v : float list
            1-D potential array
        counter : int
            returns len of v and x
        """
    x, v, counter = periodic(v_0, x_0, L, N, x_min, x_max, wells)
    eigen_simple = find_zeros(u_0, u_1, v, counter, x_min, x_max, N, e_min, e_max, accuracy)
    eigen_correct = correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy, start, eigen_simple)
    return eigen_correct, x, v, counter


def plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, graphs=0):
    
    global canvas
    
    """Calculate and plot eigenstates with given eigenenergies

            Parameters
            ----------
            x : float list
                x-axis in nm
            v : float list
                1-D potential array
            eigen_correct : float list
                list of eigenenergies in [eV]
            u_0 : float
                first boundary value
            u_1 : float
                second fixed point for numerov algorithm
            x_min : float
                start of x-axis
            x_max : float
                end of x-axis
            N : int
                steps for resolution

            Returns
            -------
            """
    fig, ax = plt.subplots()
    ax.plot(dpi=300)
    ax.plot(x, v, c="black")
    ax2 = ax.twinx()
    if graphs == 0:
        for i in range(len(eigen_correct)):
            K = get_k(v, eigen_correct[i])
            wave = numerov(u_0, u_1, K, counter, x_min, x_max)
            wave = norm(wave)
            for j, item in enumerate(wave):
                wave[j] = wave[j] + eigen_correct[i]
            ax2.plot(x, wave, label=r"$\psi_{%i}$"%(i+1))
            ax2.axhline(eigen_correct[i], ls='--', c='red')
    else:
        for i in range(graphs):
            K = get_k(v, eigen_correct[i])
            wave = numerov(u_0, u_1, K, counter, x_min, x_max)
            wave = norm(wave)
            for j, item in enumerate(wave):
                wave[j] = wave[j] + eigen_correct[i]
            ax2.plot(x, wave, label=r"$\psi_{%i}$"%(i+1))
            ax2.axhline(eigen_correct[i], ls='--', c='red')
    plt.grid(True)
    ax.set_xlabel(r'x in $[nm]$')
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red')
    plt.legend()
    
    # create a widget to show the plot in the GUI
    
    canvas = FigureCanvasTkAgg(fig, window)
    
    canvas.get_tk_widget().pack(side = tk.RIGHT, anchor="ne")
    return

def plot_bands(v_0, x_0, u_0, u_1, L, N, x_min, x_max, e_min, e_max, accuracy, max_bound, start, max_wells):
    """Calculate and plot eigenenergies for varying amount of wells

            Parameters
            ----------
            v_0 : float
                depth of potential wells in eV; inverse wells (walls) aren't supported right now
            x_0 : float
                point where the first well begins
            u_0 : float
                first boundary value
            u_1 : float
                second fixed point for numerov algorithm
            L : float
                width of the well
            N : int
                steps for resolution
            x_min : float
                start of x-axis
            x_max : float
                end of x-axis
            e_min : float
                starting point for energy search
            e_max : float
                ending point for energy search
            accuracy : float
                accuracy of search for zeros
            max_bound : int
                maximal amount of iterations if control doesn't get low enough
            start : float
                first step to calculate the new energy
            max_wells : int
                number of wells for which eigenenergies should be calculated

            Returns
            -------
            eigen_correct : float list
                list of eigenenergies in [eV]
            x : float list
                x-axis in nm
            v : float list
                1-D potential array
            counter : int
                returns len of v and x
            """
    fig, ax = plt.subplots()
    ax.plot(dpi=300)
    for i in range(1, max_wells+1):
        eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, i, e_min,
                                                             e_max, accuracy, max_bound, start)
        for j in range(len(eigen_correct)):
            ax.hlines(y=eigen_correct[j], xmin=i-1, xmax=i, linewidth=2, color='r')
    plt.grid(True)
    plt.show()


#Define the "cal" function to set the variables and start the caclulations

def cal() : 
    global on, E_max,u_0,u_1,max_bound,E_0,N,start,v_0,x_0,L,x_min,x_max,wells,E,accuracy,e_min,e_max,max_wells, on 
    global entry_Emax, entry_start, entry_N, entry_E_0, entry_u_1, entry_u_0
    if on : 
        var_v_0 = entry_tiefe.get()
        var_x_0 = entry_abstand.get()
        var_L = entry_breite.get()
        var_x_min = entry_xmin.get()
        var_x_max = entry_xmax.get()
        var_wells = entry_töpfe.get()
        var_accuracy = entry_accuracy.get()
        var_e_min = entry_emin.get()
        var_e_max = entry_emax.get()
        var_max_wells = entry_maxwells.get()
        
        u_1 = .01
        u_0 = .0
        E_max = 0
        on = 1
        E_0 = -1
        max_bound = 1000000
        start = 0.01
        N=1000
        
        
        v_0 = float(var_v_0)
        x_0 = float(var_x_0)
        L = float(var_L)
        x_min =float(var_x_min)
        x_max = float(var_x_max)
        wells = int(var_wells)
        accuracy = float(var_accuracy)
        e_min = float(var_e_min)
        e_max = float(var_e_max)
        max_wells = int(var_max_wells)
        
        eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,accuracy, max_bound, start)
        plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, 2)
        #plot_bands(v_0, x_0, u_0, u_1, L, N, x_min, x_max, e_min, e_max, accuracy, max_bound, start, max_wells)
        
        return
    else :
        var_v_0 = entry_tiefe.get()
        var_x_0 = entry_abstand.get()
        var_L = entry_breite.get()
        var_N = entry_N.get()
        var_x_min = entry_xmin.get()
        var_x_max = entry_xmax.get()
        var_wells = entry_töpfe.get()
        var_E = entry_E.get()
        var_u_0 = entry_u_0.get()
        var_u_1 = entry_u_1.get()
        var_max_bound = entry_max_bound.get()
        var_accuracy = entry_accuracy.get()
        var_start = entry_start.get()
        var_E_0 = entry_E_0.get()
        var_E_max = entry_Emax.get()
        var_e_min = entry_emin.get()
        var_e_max = entry_emax.get()
        var_max_wells = entry_maxwells.get()
    
    
        v_0 = float(var_v_0)
        x_0 = float(var_x_0)
        L = float(var_L)
        N = float(var_N)
        x_min =float(var_x_min)
        x_max = float(var_x_max)
        wells = int(var_wells)
        E = float(var_E)
        u_0 = float(var_u_0)
        u_1 = float(var_u_1)
        max_bound = int(var_max_bound)
        accuracy = float(var_accuracy)
        start = float(var_start)
        E_0 = float(var_E_0)
        E_max = float(var_E_max)
        e_min = float(var_e_min)
        e_max = float(var_e_max)
        max_wells = int(var_max_wells)
        on = 0
        eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,accuracy, max_bound, start)
        plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, 2)
        #plot_bands(v_0, x_0, u_0, u_1, L, N, x_min, x_max, e_min, e_max, accuracy, max_bound, start, max_wells)
        return
        

    
# Define the adv funtction to show and hide advance options 

def adv() : 
    global on, E_max,u_0,u_1,max_bound,E_0,N,start
    global frame_Emax, frame_start, frame_N, frame_E_0, frame_u_1, frame_u_0,frame_max_bound
    global entry_Emax, entry_start, entry_N, entry_E_0, entry_u_1, entry_u_0,entry_max_bound
    
    
    if on :
        
        frame_start = tk.Frame(master = window,bd = 5)
        label_start = tk.Label(master = frame_start, text = "start(Default is 0.01", width = 30)
        label_start.pack(side = tk.LEFT)

        entry_start = tk.Entry(master = frame_start, width = 20)
        entry_start.pack(side = tk.LEFT)
    
        frame_start.pack(anchor="nw")
        
        frame_Emax = tk.Frame(master = window,bd = 5)
        label_Emax = tk.Label(master = frame_Emax, text = "Emax(Default is 0)", width = 30)
        label_Emax.pack(side = tk.LEFT)

        entry_Emax = tk.Entry(master = frame_Emax, width = 20)
        entry_Emax.pack(side = tk.LEFT)
    
        frame_Emax.pack(anchor="nw")
        
        frame_E_0 = tk.Frame(master = window,bd = 5)
        label_E_0 = tk.Label(master = frame_E_0, text = "E_0(Default is -1)", width = 30)
        label_E_0.pack(side = tk.LEFT)

        entry_E_0 = tk.Entry(master = frame_E_0, width = 20)
        entry_E_0.pack(side = tk.LEFT)
    
        frame_E_0.pack(anchor="nw")
        
        frame_N = tk.Frame(master = window,bd = 5)
        label_N = tk.Label(master = frame_N, text = "N(Default is 1000 )", width = 30)
        label_N.pack(side = tk.LEFT)

        entry_N = tk.Entry(master = frame_N, width = 20)
        entry_N.pack(side = tk.LEFT)
    
        frame_N.pack(anchor="nw")
        
        frame_max_bound = tk.Frame(master = window,bd = 5)
        label_max_bound = tk.Label(master = frame_max_bound, text = "max_bound(Default is 1000000 )", width = 30)
        label_max_bound.pack(side = tk.LEFT)

        entry_max_bound = tk.Entry(master = frame_max_bound, width = 20)
        entry_max_bound.pack(side = tk.LEFT)
    
        frame_max_bound.pack(anchor="nw")
        
        frame_u_0 = tk.Frame(master = window,bd = 5)
        label_u_0 = tk.Label(master = frame_u_0, text = "u_0(Default is .0", width = 30)
        label_u_0.pack(side = tk.LEFT)

        entry_u_0 = tk.Entry(master = frame_u_0, width = 20)
        entry_u_0.pack(side = tk.LEFT)
    
        frame_u_0.pack(anchor="nw")
        
        frame_u_1 = tk.Frame(master = window,bd = 5)
        label_u_1 = tk.Label(master = frame_u_1, text = "u_1(Defualt is .01", width = 30)
        label_u_1.pack(side = tk.LEFT)

        entry_u_1 = tk.Entry(master = frame_u_1, width = 20)
        entry_u_1.pack(side = tk.LEFT)
    
        frame_u_1.pack(anchor="nw")
        on = 0
    else :
        frame_E_0.destroy()
        frame_N.destroy()
        frame_u_1.destroy()
        frame_Emax.destroy()
        frame_u_0.destroy()
        frame_start.destroy()
        frame_max_bound.destroy()
        u_1 = .01
        u_0 = .0
        E_max = 0
        on = 1
        E_0 = -1
        max_bound = 1000000
        start = 0.01
        N=1000

def clear() : 
    
    canvas.get_tk_widget().destroy()
    return

# set the size of the window 


window.geometry("1200x700")

#create Title of the window

window.title("Periodische Abfolge endlicher Potentialtöpfe")

# create frames for each variable in the project

frame_Graph = tk.Frame()


frame_töpfe = tk.Frame(master = window, bd= 5)

frame_tiefe = tk.Frame(master = window, bd = 5)

frame_breite = tk.Frame(master = window,bd = 5)

frame_abstand = tk.Frame(master = window,bd = 5)

frame_knöpfe = tk.Frame(master = window,bd = 10)

frame_emin = tk.Frame(master = window,bd = 5)

frame_emax = tk.Frame(master = window,bd = 5)

frame_xmin = tk.Frame(master = window,bd = 5)

frame_xmax = tk.Frame(master = window,bd = 5)

frame_E = tk.Frame(master = window,bd = 5)

frame_maxwells = tk.Frame(master = window,bd = 5)

frame_accuracy = tk.Frame(master = window,bd = 5)

frame_advance = tk.Frame(master = window,bd = 5)

# create the label for the number of the pots with an entry

label_töpfe = tk.Label(master = frame_töpfe, text = "Number of wells", width = 30)
label_töpfe.pack(side = tk.LEFT)

entry_töpfe = tk.Entry(master = frame_töpfe, width = 20)
entry_töpfe.pack(side = tk.LEFT)

# create the label for the wight of the pots with an entry

label_breite = tk.Label(master = frame_breite, text = "Width of the wells(in nm)", width = 30)
label_breite.pack(side = tk.LEFT)


entry_breite = tk.Entry(master = frame_breite, width = 20)
entry_breite.pack(side = tk.LEFT)

# create the label for the depth of the pots  with an entry

label_tiefe = tk.Label(master = frame_tiefe, text = "Depth of the wells(in nm)", width = 30)
label_tiefe.pack(side = tk.LEFT)

entry_tiefe = tk.Entry(master = frame_tiefe, width = 20)
entry_tiefe.pack(side = tk.LEFT)

# create a button to clalculate and close the window
# and a button to clear the graphs 

button_calc = ttk.Button(master = frame_knöpfe,text = "Exit", command= window.destroy)
button_calc.pack(side= tk.RIGHT,)

button_calc = ttk.Button(master = frame_knöpfe,text = "Calculate", command= cal)
button_calc.pack(side= tk.RIGHT,)

button_clear = ttk.Button(master = frame_knöpfe,text = "clear", command= clear)
button_clear.pack(side= tk.RIGHT,)


#create the label for distance between the pots 

label_abstand = tk.Label(master = frame_abstand, text = "Distance between the wells(in nm)", width = 30)
label_abstand.pack(side = tk.LEFT)

entry_abstand = tk.Entry(master = frame_abstand, width = 20)
entry_abstand.pack(side = tk.LEFT)

#create the label for emin

label_emin = tk.Label(master = frame_emin, text = "emin", width = 30)
label_emin.pack(side = tk.LEFT)

entry_emin = tk.Entry(master = frame_emin, width = 20)
entry_emin.pack(side = tk.LEFT)

# create the label for emax

label_emax = tk.Label(master = frame_emax, text = "emax", width = 30)
label_emax.pack(side = tk.LEFT)

entry_emax = tk.Entry(master = frame_emax, width = 20)
entry_emax.pack(side = tk.LEFT)

# create the label for E

label_E = tk.Label(master = frame_E, text = "E", width = 30)
label_E.pack(side = tk.LEFT)

entry_E = tk.Entry(master = frame_E, width = 20)
entry_E.pack(side = tk.LEFT)

# create the label for accuracy

label_accuracy = tk.Label(master = frame_accuracy, text = "accuracy", width = 30)
label_accuracy.pack(side = tk.LEFT)

entry_accuracy = tk.Entry(master = frame_accuracy, width = 20)
entry_accuracy.pack(side = tk.LEFT)

# create the label for maxwells

label_maxwells = tk.Label(master = frame_maxwells, text = "max wells", width = 30)
label_maxwells.pack(side = tk.LEFT)

entry_maxwells = tk.Entry(master = frame_maxwells, width = 20)
entry_maxwells.pack(side = tk.LEFT)

# create the label for xmin

label_xmin = tk.Label(master = frame_xmin, text = "xmin", width = 30)
label_xmin.pack(side = tk.LEFT)

entry_xmin = tk.Entry(master = frame_xmin, width = 20)
entry_xmin.pack(side = tk.LEFT)

# create the label for xmax

label_xmax = tk.Label(master = frame_xmax, text = "xmax", width = 30)
label_xmax.pack(side = tk.LEFT)

entry_xmax = tk.Entry(master = frame_xmax, width = 20)
entry_xmax.pack(side = tk.LEFT)

# create checkbox for advance GUI
Button_advance = ttk.Button(master = frame_advance, text = "Advance Options", command = adv)                               
Button_advance.pack(side = tk.LEFT)

#pack all the frames, buttons and labels

frame_töpfe.pack( anchor="nw" )

frame_tiefe.pack(anchor="nw") 

frame_breite.pack(anchor="nw")

frame_abstand.pack(anchor="nw")

frame_accuracy.pack(anchor="nw")

frame_emin.pack(anchor="nw")

frame_emax.pack(anchor="nw")

frame_xmin.pack(anchor="nw")

frame_xmax.pack(anchor="nw")

frame_E.pack(anchor="nw")

frame_maxwells.pack(anchor="nw")

frame_advance.pack(anchor="nw")

frame_Graph.pack(anchor="ne")

frame_knöpfe.pack(side = tk.BOTTOM, anchor="se")




window.mainloop()

var_v_0 = entry_tiefe.get()
var_x_0 = entry_abstand.get()
var_L = entry_breite.get()
var_N = entry_N.get()
var_x_min = entry_xmin.get()
var_x_max = entry_xmax.get()
var_wells = entry_töpfe.get()
var_E = entry_E.get()
var_u_0 = entry_u_0.get()
var_u_1 = entry_u_1.get()
var_max_bound = entry_max_bound.get()
var_accuracy = entry_accuracy.get()
var_start = entry_start.get()
var_E_0 = entry_E_0.get()
var_E_max = entry_Emax.get()
var_e_min = entry_emin.get()
var_e_max = entry_emax.get()
var_max_wells = entry_maxwells.get()


v_0 = int(var_v_0)
x_0 = int(var_x_0)
L = int(var_L)
N = float(var_N)
x_min =int(var_x_min)
x_max = int(var_x_max)
wells = int(var_wells)
E = int(var_E)
u_0 = int(var_u_0)
u_1 = int(var_u_1)
max_bound = int(var_max_bound)
accuracy = int(var_accuracy)
start = int(var_start)
E_0 = int(var_E_0)
E_max = int(var_E_max)
e_min = int(var_e_min)
e_max = int(var_e_max)
max_wells = int(var_max_wells)

# start the mainloop







