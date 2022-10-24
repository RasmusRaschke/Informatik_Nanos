# Sample Code zur Vorlesung Informatik für Nanowissenschaftler an der Universität Hamburg WS 20/21
# Autor: Prof. Dr. Rübhausen
# Dieser Code dient zu Illustration der Vorlesung. Schreiben Sie ihren eigenen Code für ihre Projekte - bitte kein "Copy & Paste" sie reduzieren
# sonst ihren eigenen Lerneffekt und schaden sich selbst


import matplotlib.pyplot as plt # Erzeugt eine Instanz von Pyplot
import numpy as np # Erzeugt eine Instanz von Numpy

# Plotten einer Referenzfunktion
def sin(K,Delta_X,N):# Definiere eine Prozedur um mit den sin darzustellen
    Werte_Sin = []#Liste mit Werten des Sinus
    h = Delta_X/N # Schrittweite
    print("Schrittweite")#Ausgabe
    print(h)
    for i in range(1,N):#Schleife zur Berechnung
        Werte_Sin.append(np.sin(K*i*h))#Numpy Sin
        
    return Werte_Sin #Rückgabe der Lösung



# numerische Lösungs-Prozeduren
def Euler_DGL(K,U0,Z0,Delta_X,N):# Definiere eine Prozedur um mit dem Euler-Verfahren eine Differentialgleichung 2ter Ordnung auszurechnen. Die Parameter der Prozedur ... K = omega^2
    Werte_U = []#Liste mit Werten der Lösung der DGL
    Werte_Z = []# List mit Werten der Ableitung der Lösung z.b. bei einem Federpendel Ort und Geschwindigkeit
    Werte_U.append(U0) # Trägt den Start von U in den Listenplatz "Null" ein
    Werte_Z.append(Z0) # Trägt den Start von Z in den Listenplatz "Null ein
    h = Delta_X/N # Schrittweite
    print("Schrittweite")#Ausgabe
    print(h)
    for i in range(1,N):#Schleife zur Berechnung
        Werte_U.append(Werte_Z[i-1]*h+Werte_U[i-1])#Euler DGL1 für DGL 1 Ordnung
        Werte_Z.append((-K)*Werte_U[i-1]*h+Werte_Z[i-1])#Euler DGL 2 für DGL 1 Ordnung
    return Werte_U#Rückgabe der Lösung

def Num_DGL(K,U0,U1,Delta_X,N):# Definiere eine Prozedur um mit dem Numerov -Verfahren eine Differentialgleichung 2ter Ordnung auszurechnen. Die Parameter der Prozedur ... K = omega^2
    Werte_U = []#s.o.
    Werte_U.append(U0) # Trägt den Start von U in den Listenplatz "Null" ein
    Werte_U.append(U1) # Wir benötign den Wert x-h und x um x+h abzuschätzen
    h = Delta_X/N
    print("Schrittweite")
    print(h)
    for i in range(2,N):
        Werte_U.append(((2 * Werte_U[i - 1] * (1 - 5 / 12 * h * h * K)) - Werte_U[i - 2] * (
                    1 + h * h * K / 12)) / (1 + 1 / 12 * h * h * K))  # Numerov Approximation mit konstantem K
    return Werte_U



# Parameter - global

Steps = 1000 # Schritte
Werte_Bereich = 100 # Zeitraum oder Distanz
Start_Gesch = 1
Start_Auslenkung = 0
Frequenz_zumQuadrat = 1
Num_2 = 0
Num_1 = 0.1

Werte_Sin = sin(np.sqrt(Frequenz_zumQuadrat),Werte_Bereich,Steps) # hier muss ja sin(omega t) stehen nicht omega^2t

Lösung_Euler = Euler_DGL(Frequenz_zumQuadrat,Start_Auslenkung,Start_Gesch,Werte_Bereich,Steps)
Lösung_Numerov = Num_DGL(Frequenz_zumQuadrat,Num_2,Num_1,Werte_Bereich,Steps)


plt.plot(Werte_Sin) # Man sieht 2pi/omega !
#plt.plot(Lösung_Euler)
plt.plot(Lösung_Numerov)



plt.show()
