#!/usr/bin/env python3

import matplotlib.pyplot as plt # Erzeugt eine Instanz von Pyplot
import numpy as np # Erzeugt eine Instanz von Numpy

# Definition einer Potentialwandlandschaft - Tiefe V_0, (eV) Position x_0 (nm) und Breite L (nm) - Delta_X ist der Bereich indem gerechnet wird (nm), N = Anzahl der Schritte in x

def Potential_Landschaft(V0,X0,L,Delta_X,N):
	Werte_Potential = []
	X_Achse = []
	h = Delta_X / N # Auflösung in nm /schritt
	print("Auflösung in nm")
	print(h)
	Schritte_bis_X0 = round(X0/h)#Ganze Zahl für Schleife
	print("Schritte bis X0")
	print(Schritte_bis_X0)
	Schritte_ohne_Potential = round(L/h)#Ganze Zahl für Schleife
	print("Schritte zwischen den Wänden")
	print(Schritte_ohne_Potential)
	Schritte_nach_Potential = N-Schritte_bis_X0 - Schritte_ohne_Potential # Brauchen wir nicht
	
	for i in range(N):
		X_Achse.append(i*h)#Erzeugung der x-Achse in nm
		
	for i in range(Schritte_bis_X0):#Bis X0 gibt es eine Potentialwand
		Werte_Potential.append(V0)
	for i in range(Schritte_bis_X0,(Schritte_bis_X0 + Schritte_ohne_Potential)):#Der Bereich ab X0 im Bereich L ist potentialfrei
		Werte_Potential.append(0)
	for i in range((Schritte_bis_X0 + Schritte_ohne_Potential),N):#Danach gibt es die rechte Potentialwand
		Werte_Potential.append(V0)
	
	return Werte_Potential, X_Achse#Rückgabe des Potentials (eV) und der X-Achse (nm)

# jetzt müssen wir noch die richtige K_List berechnen - dass ist der Omega^2 Vorfaktor und er hängt von der Energie des Teilchens ab

def K_list_gen(Potential,Energie,N):
	K_List = []
	Faktor = 26.27#Umrechnungsfaktor für eV und nm
	for i in range(N):
		K_List.append(Faktor * (Energie-Potential[i]))#Die Energie des Teilchens ist fest aber die Potential-Landschaft ändert sich mit der Position in X
	return K_List # Rückgabe des ortsabhängigen K

def Wavefunction(K,U0,U1,Delta_X,N):#Berechnung der Wellenfunktion mit ortsabhängigem K
	Werte_U = []
	Werte_U.append(U0)  # Trägt den Start von U in den Listenplatz "Null" ein
	Werte_U.append(U1)  # Wir benötign den Wert x-h und x um x+h abzuschätzen
	h = Delta_X / N
	# print("Schrittweite")
	# print(h)
	for i in range(2, N):
		Werte_U.append(((2 * Werte_U[i - 1] * (1 - 5 / 12 * h * h * K[i-1])) - Werte_U[i - 2] * (
				1 + h * h * K[i-2] / 12)) / (1 + 1 / 12 * h * h * K[i]))  # Numerov mit K(x) vergleiche mit Numerov oben
		
	return Werte_U

def Berechnung_der_Normierten_Wellenfunktion(Eigenwert,U0,U1,Potential,Delta_X,N):
	Norm_Wave =[]
	frequency = K_list_gen(Potential,Eigenwert,N)
	Eigen_Wave = Wavefunction(frequency,U0,U1,Delta_X,N)
	# Wir haben jetzt eine nicht-normierte Wellenfunktion zum Energie-Eigenwert ausgerechnet
	# Wir müssen jetzt die Normierung ausrechnen
	Summe = 0
	
	for i in range(N):#Summe = Integral
		Summe = Eigen_Wave[i]**2+Summe
	print("Normierung",Summe)
	for i in range(N):
		Norm_Wave.append((Eigen_Wave[i]/np.sqrt(Summe))) #Normierung der Werte
		
	return Norm_Wave


Num_0 = 0
Num_1 = 0.000001
Range = 15 # in x von 0 bis 15 nm
Position = 2.5 # Position linke Wand
Length =  10 # Ohne Potential
Steps = 1000 # Schritte
Potential_Wert = 1 #  Potential in eV

# 1 eV (compare to 0.00375637 (inf) / 0.0034764667 (n=1) / 0.0139058668 n=2 (nom.) 0.013904041

Pot, X =Potential_Landschaft(Potential_Wert,Position,Length,Range,Steps)

Wave_1 = Berechnung_der_Normierten_Wellenfunktion(0.0034764667,Num_0,Num_1,Pot,Range,Steps)
Wave_2 = Berechnung_der_Normierten_Wellenfunktion(0.013904041005392142,Num_0,Num_1,Pot,Range,Steps)
Wave_3 = Berechnung_der_Normierten_Wellenfunktion(0.03127716244972,Num_0,Num_1,Pot,Range,Steps)

fig, ax1 = plt.subplots()

color = "tab:red"
ax1.set_xlabel("X-Axis (nm)")
ax1.set_ylabel("Potential (eV)")
ax1.plot(X, Pot, color = color)
ax1.tick_params(axis='y', labelcolor = color)


ax2 = ax1.twinx()

color = "tab:blue"
ax2.set_ylabel("Psi", color = color)
#ax2.plot(X, Wave, color = color)
ax2.plot(X, Wave_1, color = color)
#ax2.plot(X, Wave_2, color=color)
#ax2.plot(X, Wave_3, color = color)
ax2.tick_params(axis='y', labelcolor =color)

fig.tight_layout()

#plt.xlim(1,14) #Achsen einschränken
#plt.ylim(-50,50)

plt.show()