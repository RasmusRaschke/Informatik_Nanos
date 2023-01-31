import numpy as np
import matplotlib.pyplot as plt
import math as m
plt.rc('legend', fontsize=8)
########################################################################################################################
def periodic(v_0, x_0, L, N, x_min, x_max, charge, el_field, wells=1, electric=1):
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
    if electric == 1:
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
    else:
        for i in range(wells):
            for j in range(spacing_steps):
                v.append(charge * el_field * j)
                counter += 1
            for k in range(well_steps):
                v.append(- v_0 + charge * el_field * k)
                counter += 1
        for i in range(spacing_steps):
            v.append(charge * el_field * i)
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
    print("e_old:", e)
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
        control = abs(e_new - e)
        print("c:", control)
        e = e_new
        if control < accuracy:
            break
    print("e_new:", e_new)
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


def calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max, accuracy_cheap, accuracy_exp,
                          max_bound, start, charge, el_field, electric=1):
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
    x, v, counter = periodic(v_0, x_0, L, N, x_min, x_max, charge, el_field, wells, electric)
    eigen_simple = find_zeros(u_0, u_1, v, counter, x_min, x_max, N, e_min, e_max, accuracy_cheap)
    eigen_correct = correct_zeros(u_0, u_1, counter, x_min, x_max, N, v, max_bound, accuracy_exp, start, eigen_simple)
    return eigen_correct, x, v, counter


def plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, filename, legend=0, graphs=0):
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
            graphs : int
                decide if all eigenfunctions or only certain amount should be plotted

            Returns
            -------
            """
    fig, ax = plt.subplots()
    ax.plot(dpi=1000)
    ax.plot(x, v, c="black")
    ax2 = ax.twinx()
    if graphs == 0:
        for i in range(len(eigen_correct)):
            K = get_k(v, eigen_correct[i])
            wave = numerov(u_0, u_1, K, counter, x_min, x_max)
            wave = norm(wave)
            for j, item in enumerate(wave):
                wave[j] = wave[j] + eigen_correct[i]
            ax2.plot(x, wave, label=r"$\varepsilon_{%i} = %f $ eV"%(i+1, eigen_correct[i]))
            ax2.axhline(eigen_correct[i], ls='--', c='red')
    else:
        for i in range(graphs):
            K = get_k(v, eigen_correct[i])
            wave = numerov(u_0, u_1, K, counter, x_min, x_max)
            wave = norm(wave)
            for j, item in enumerate(wave):
                wave[j] = wave[j] + eigen_correct[i]
            ax2.plot(x, wave, label=r"$\varepsilon_{%i} = %f $ eV"%(i+1, eigen_correct[i]))
            ax2.axhline(eigen_correct[i], ls='--', c='red')
    plt.grid(True)
    ax.set_xlabel(r'$x$ in $[nm]$')
    ax.set_ylabel(r'$v$ in $[eV]$')
    ax2.set_ylabel(r'$\Psi$', c='red')
    ax2.spines['right'].set_color('red')
    ax2.tick_params(axis='y', colors='red')
    plt.legend(loc='upper left')
    if legend == 1:
        plt.legend('', frameon=False)
    plt.savefig('/home/rasmus/Informatik_Nanos/plots/%s.png'%filename, dpi=1000)
    plt.show()


def plot_bands(v_0, x_0, u_0, u_1, L, N, x_min, x_max, e_min, e_max, accuracy_cheap, accuracy_exp, max_bound,
               start, max_wells, charge, el_field, electric, filename):
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
            accuracy_cheap : float
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
    ax.plot(dpi=1200)
    for i in range(1, max_wells+1):
        eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, i, e_min, e_max,
                                                    accuracy_cheap, accuracy_exp, max_bound, start, charge, el_field,
                                                     electric)
        for j in range(len(eigen_correct)):
            ax.hlines(y=eigen_correct[j], xmin=i-1, xmax=i, linewidth=2, color='r')
    plt.grid(True)
    ax.set_xlabel(r'Kastenanzahl')
    ax.set_ylabel(r'Eigenenergien in $[eV]$')
    plt.savefig('/home/rasmus/Informatik_Nanos/plots/%s.png'%filename, dpi=1000)
    plt.show()


########################################################################################################################
#Variablen hier
v_0 = 8 #Kastentiefe in [eV]
x_0 = 2 #Startpunkt des ersten Kastens in [nm]
L = 0.287 #Kastenbreite in [nm]; Silber: 0.287
wells = 2 #Kastenanzahl
accuracy_cheap = 0.01 #Genauigkeit
e_min = -8 #untere Grenze der Eigenwertsuche
e_max = 15 #obere Grenze der Eigenwertsuche
max_wells = 200 #Kastenanzahl für Bandberechnung
filename = 'test' #NICHT IN GUI
legend = 0 #NICHT IN GUI
charge = 0.30282212 #Ladung des Teilchens [einheitenlos]
el_field = 0.05 #Elektrisches Feld in [eV^2]
electric = 1 #DAS SOLL EIN WAHLFELD WERDEN, ALSO ELEKTISCHES FELD JA/NEIN ANKREUZEN
E = 1. #ENTWICKLEROPTION; Energie der k-List in [eV]
u_0 = .0 #ENTWICKLEROPTION; 1. Randbedingung des Numerov-Verfahrens
u_1 = .001 #ENTWICKLEROPTION; 2. Randbedingung des Numerov-Verfahrens
max_bound = 1000000 #ENTWICKLEROPTION; maximale Iterationen des NR-Verfahrens
accuracy_exp = 0.00000001 #ENTWICKLEROPTION; Genauigkeit des NR-Verfahrens
start = 0.000000001 #ENTWICKLEROPTION; erstes Intervall der Eigenwertsuche
N = 1000 #ENTWICKLEROPTION #Auflösung des Potentials
x_min = 0 #ENTWICKLEROPTION; Beginn des Plots
x_max = 10 #ENTWICKLEROPTION; Ende des Plots
########################################################################################################################
#Berechnung #SO FUNKTIONIERT DAS PROGRAMM MIT DEN AKTUELLEN VARIABLEN
eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,
                                                    accuracy_cheap, accuracy_exp, max_bound, start, charge, el_field,
                                                    electric)
#plot_eigenstates(x, v, eigen_correct, u_0, u_1, counter, x_min, x_max, filename, legend)

plot_bands(v_0, x_0, u_0, u_1, L, N, x_min, x_max, e_min, e_max, accuracy_cheap, accuracy_exp, max_bound, start,
           max_wells, charge, el_field, electric, filename)

########################################################################################################################
#EINFACH IGNORIEREN
'''
average_x = [1, 1.5, 2, 2.5, 3, 5, 10]
plot_average_1 = [-0.2559052422454842, -0.34271239128133507, -0.38613615166845133, -0.4105402983419828, 
-0.4257988851412097, -0.45321103550772684, -0.47045380448876606]
plot_average_2 = [-0.29137048925821984, -0.3546134550060287, -0.38809011653195064, -0.40789164055927113, 
-0.4207453397473255, -0.4450937912244151, -0.4617850941578363]
plot_average_3 = [-0.17527399182584188, -0.22804813899323487, -0.26463821991839376, -0.29132183753909535, 
-0.3116466744996986, -0.36061326839469215, -0.4090776227071717]
plt.scatter(average_x, plot_average_1)
plt.scatter(average_x, plot_average_2)
plt.scatter(average_x, plot_average_3)
plt.plot(average_x, plot_average_1, label="Ein Kasten")
plt.plot(average_x, plot_average_2, label="Zwei Kästen")
plt.plot(average_x, plot_average_3, label="Zehn Kästen")
plt.legend(loc='upper right')
plt.grid(True)
plt.xlabel(r"Topfabstand in $[nm]$")
plt.ylabel(r'Mittel der ersten drei Eigenenergien in $[eV]$')
plt.savefig('/home/rasmus/Informatik_Nanos/plots/width_averages.png', dpi=1200)
plt.show()
'''

'''
fign, axn = plt.subplots()
axn.plot(dpi=1200)
plt.grid(True)
axn.set_xlabel(r'elektrische Feldstärke in $[eV^2]$')
axn.set_ylabel(r'Eigenenergien in $[eV]$')
average1 = []
average2 = []
average3 = []
field = [0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01]
for j in field:
    v_0 = 5  # Kastentiefe in [eV]
    x_0 = 2  # Startpunkt des ersten Kastens in [nm]
    L = 0.287  # Kastenbreite in [nm]
    wells = 1  # Kastenanzahl
    accuracy_cheap = 0.001  # Genauigkeit
    e_min = -5  # untere Grenze der Eigenwertsuche
    e_max = 0  # obere Grenze der Eigenwertsuche
    max_wells = 10  # Kastenanzahl für Bandberechnung
    filename = 'electric_field_%f'%j  # NICHT IN GUI
    legend = 0  # NICHT IN GUI
    charge = 0.30282212  # Ladung des Teilchens [einheitenlos]
    el_field = j  # Elektrisches Feld in [eV^2]
    electric = 0  # DAS SOLL EIN WAHLFELD WERDEN, ALSO ELEKTISCHES FELD JA/NEIN ANKREUZEN
    E = 1.  # ENTWICKLEROPTION; Energie der k-List in [eV]
    u_0 = .0  # ENTWICKLEROPTION; 1. Randbedingung des Numerov-Verfahrens
    u_1 = .001  # ENTWICKLEROPTION; 2. Randbedingung des Numerov-Verfahrens
    max_bound = 1000000  # ENTWICKLEROPTION; maximale Iterationen des NR-Verfahrens
    accuracy_exp = 0.00000001  # ENTWICKLEROPTION; Genauigkeit des NR-Verfahrens
    start = 0.000000001  # ENTWICKLEROPTION; erstes Intervall der Eigenwertsuche
    N = 1000  # ENTWICKLEROPTION #Auflösung des Potentials
    x_min = 0  # ENTWICKLEROPTION; Beginn des Plots
    x_max = 10  # ENTWICKLEROPTION; Ende des Plots
    eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,
                                                         accuracy_cheap, accuracy_exp, max_bound, start, charge,
                                                         el_field,
                                                         electric)
    average = 0
    for k in range(2):
        average += eigen_correct[k]
    average = average / 2
    average1.append(average)

for j in field:
    v_0 = 5  # Kastentiefe in [eV]
    x_0 = 2  # Startpunkt des ersten Kastens in [nm]
    L = 0.287  # Kastenbreite in [nm]
    wells = 3  # Kastenanzahl
    accuracy_cheap = 0.001  # Genauigkeit
    e_min = -5  # untere Grenze der Eigenwertsuche
    e_max = 0  # obere Grenze der Eigenwertsuche
    max_wells = 10  # Kastenanzahl für Bandberechnung
    filename = 'electric_field_%f'%j  # NICHT IN GUI
    legend = 0  # NICHT IN GUI
    charge = 0.30282212  # Ladung des Teilchens [einheitenlos]
    el_field = j  # Elektrisches Feld in [eV^2]
    electric = 0  # DAS SOLL EIN WAHLFELD WERDEN, ALSO ELEKTISCHES FELD JA/NEIN ANKREUZEN
    E = 1.  # ENTWICKLEROPTION; Energie der k-List in [eV]
    u_0 = .0  # ENTWICKLEROPTION; 1. Randbedingung des Numerov-Verfahrens
    u_1 = .001  # ENTWICKLEROPTION; 2. Randbedingung des Numerov-Verfahrens
    max_bound = 1000000  # ENTWICKLEROPTION; maximale Iterationen des NR-Verfahrens
    accuracy_exp = 0.00000001  # ENTWICKLEROPTION; Genauigkeit des NR-Verfahrens
    start = 0.000000001  # ENTWICKLEROPTION; erstes Intervall der Eigenwertsuche
    N = 1000  # ENTWICKLEROPTION #Auflösung des Potentials
    x_min = 0  # ENTWICKLEROPTION; Beginn des Plots
    x_max = 10  # ENTWICKLEROPTION; Ende des Plots
    eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,
                                                         accuracy_cheap, accuracy_exp, max_bound, start, charge,
                                                         el_field,
                                                         electric)
    average = 0
    for k in range(2):
        average += eigen_correct[k]
    average = average / 2
    average2.append(average)

for j in field:
    v_0 = 5  # Kastentiefe in [eV]
    x_0 = 2  # Startpunkt des ersten Kastens in [nm]
    L = 0.287  # Kastenbreite in [nm]
    wells = 10  # Kastenanzahl
    accuracy_cheap = 0.001  # Genauigkeit
    e_min = -5  # untere Grenze der Eigenwertsuche
    e_max = 0  # obere Grenze der Eigenwertsuche
    max_wells = 10  # Kastenanzahl für Bandberechnung
    filename = 'electric_field_%f'%j  # NICHT IN GUI
    legend = 0  # NICHT IN GUI
    charge = 0.30282212  # Ladung des Teilchens [einheitenlos]
    el_field = j  # Elektrisches Feld in [eV^2]
    electric = 0  # DAS SOLL EIN WAHLFELD WERDEN, ALSO ELEKTISCHES FELD JA/NEIN ANKREUZEN
    E = 1.  # ENTWICKLEROPTION; Energie der k-List in [eV]
    u_0 = .0  # ENTWICKLEROPTION; 1. Randbedingung des Numerov-Verfahrens
    u_1 = .001  # ENTWICKLEROPTION; 2. Randbedingung des Numerov-Verfahrens
    max_bound = 1000000  # ENTWICKLEROPTION; maximale Iterationen des NR-Verfahrens
    accuracy_exp = 0.00000001  # ENTWICKLEROPTION; Genauigkeit des NR-Verfahrens
    start = 0.000000001  # ENTWICKLEROPTION; erstes Intervall der Eigenwertsuche
    N = 1000  # ENTWICKLEROPTION #Auflösung des Potentials
    x_min = 0  # ENTWICKLEROPTION; Beginn des Plots
    x_max = 10  # ENTWICKLEROPTION; Ende des Plots
    eigen_correct, x, v, counter = calculate_eigenvalues(v_0, x_0, u_0, u_1, L, N, x_min, x_max, wells, e_min, e_max,
                                                         accuracy_cheap, accuracy_exp, max_bound, start, charge,
                                                         el_field,
                                                         electric)
    average = 0
    for k in range(2):
        average += eigen_correct[k]
    average = average / 2
    average3.append(average)
axn.scatter(field, average1, color='black')
axn.plot(field, average1, label='ein Kasten', color='black')
axn.scatter(field, average2, color='red')
axn.plot(field, average2, label='drei Kästen', color='red')
axn.scatter(field, average3, color='blue')
axn.plot(field, average3, label='zehn Kästen', color='blue')
plt.legend(loc='upper left')
plt.savefig('/home/rasmus/Informatik_Nanos/plots/electric_field_strength.png', dpi=1000)
plt.show()
'''
########################################################################################################################