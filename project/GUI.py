<<<<<<< HEAD

# import modules
import tkinter as tk
from tkinter import ttk

#Define the "cal" function to chek if the entry is a number else print it to the user
def cal() : 
    if entry_töpfe.get().isdigit():
        labelr["text"] = ""
    else :
        labelr["text"] = "{Number of pots is not a Number}"
        return
    if entry_tiefe.get().isdigit():
        labelr["text"] = ""
    else :
        labelr["text"] = "{Depth is not a Number}"
        return
    if entry_breite.get().isdigit():
         labelr["text"] = ""
    else :
         labelr["text"] = "{Width is not a Number}"
         return
    if entry_abstand.get().isdigit():
          labelr["text"] = ""
    else :
          labelr["text"] = "{Distance is not a Number}"
          return                        
    return

# create window a window where everything will be shown 
window = tk.Tk()
window.geometry("600x300")
#create Title of the window
window.title("Periodische Abfolge endlicher Potentialtöpfe")

# create frames for each variable in the project
frame_töpfe = tk.Frame(bd = 5)

frame_tiefe = tk.Frame(bd = 5)

frame_breite = tk.Frame(bd = 5)

frame_abstand = tk.Frame(bd = 5)

frame_knöpfe = tk.Frame(bd = 10)

# create the label for the number of the pots with an entry
label_töpfe = tk.Label(master = frame_töpfe, text = "Number of pots", width = 20)
label_töpfe.pack(side = tk.LEFT)

entry_töpfe = tk.Entry(master = frame_töpfe, width = 20)
entry_töpfe.pack(side = tk.LEFT)

# create the label for the wight of the pots with an entry
label_breite = tk.Label(master = frame_breite, text = "Width of the pots", width = 20)
label_breite.pack(side = tk.LEFT)


entry_breite = tk.Entry(master = frame_breite, width = 20)
entry_breite.pack(side = tk.LEFT)

# create the label for the depth of the pots  with an entry
label_tiefe = tk.Label(master = frame_tiefe, text = "Depth of the plots", width = 20)
label_tiefe.pack(side = tk.LEFT)

entry_tiefe = tk.Entry(master = frame_tiefe, width = 20)
entry_tiefe.pack(side = tk.LEFT)

# create a button to clalculate and close he window

button_calc = ttk.Button(master = frame_knöpfe,text = "Exit", command= window.destroy)
button_calc.pack(side= tk.RIGHT,)

button_calc = ttk.Button(master = frame_knöpfe,text = "Calculate", command= cal)
button_calc.pack(side= tk.RIGHT,)

#create a label for the result or show when the input is wrong
labelr = tk.Label(text = "" , bd = 10)


#create the label for distance between the pots 
label_abstand = tk.Label(master = frame_abstand, text = "Distance between the pots", width = 20)
label_abstand.pack(side = tk.LEFT)

entry_abstand = tk.Entry(master = frame_abstand, width = 20)
entry_abstand.pack(side = tk.LEFT)

#pack all the frames, buttons and labels
frame_töpfe.pack( anchor="nw" )

frame_tiefe.pack(anchor="nw") 

frame_breite.pack(anchor="nw")

frame_abstand.pack(anchor="nw")

frame_knöpfe.pack(side = tk.BOTTOM, anchor="se")

labelr.pack(anchor = "nw")
    
# start the mainloop
window.mainloop()


=======

# import modules
import tkinter as tk
from tkinter import ttk

#Define the "cal" function to chek if the entry is a number else print it to the user
def cal() : 
    if entry_töpfe.get().isdigit():
        labelr["text"] = ""
    else :
        labelr["text"] = "{Number of pots is not a Number}"
        return
    if entry_tiefe.get().isdigit():
        labelr["text"] = ""
    else :
        labelr["text"] = "{Depth is not a Number}"
        return
    if entry_breite.get().isdigit():
         labelr["text"] = ""
    else :
         labelr["text"] = "{Width is not a Number}"
         return
    if entry_abstand.get().isdigit():
          labelr["text"] = ""
    else :
          labelr["text"] = "{Distance is not a Number}"
          return                        
    return

# create window a window where everything will be shown 
window = tk.Tk()
window.geometry("600x300")
#create Title of the window
window.title("Periodische Abfolge endlicher Potentialtöpfe")

# create frames for each variable in the project
frame_töpfe = tk.Frame(bd = 5)

frame_tiefe = tk.Frame(bd = 5)

frame_breite = tk.Frame(bd = 5)

frame_abstand = tk.Frame(bd = 5)

frame_knöpfe = tk.Frame(bd = 10)

# create the label for the number of the pots with an entry
label_töpfe = tk.Label(master = frame_töpfe, text = "Number of pots", width = 20)
label_töpfe.pack(side = tk.LEFT)

entry_töpfe = tk.Entry(master = frame_töpfe, width = 20)
entry_töpfe.pack(side = tk.LEFT)

# create the label for the wight of the pots with an entry
label_breite = tk.Label(master = frame_breite, text = "Width of the pots", width = 20)
label_breite.pack(side = tk.LEFT)


entry_breite = tk.Entry(master = frame_breite, width = 20)
entry_breite.pack(side = tk.LEFT)

# create the label for the depth of the pots  with an entry
label_tiefe = tk.Label(master = frame_tiefe, text = "Depth of the plots", width = 20)
label_tiefe.pack(side = tk.LEFT)

entry_tiefe = tk.Entry(master = frame_tiefe, width = 20)
entry_tiefe.pack(side = tk.LEFT)

# create a button to clalculate and close he window

button_calc = ttk.Button(master = frame_knöpfe,text = "Exit", command= window.destroy)
button_calc.pack(side= tk.RIGHT,)

button_calc = ttk.Button(master = frame_knöpfe,text = "Calculate", command= cal)
button_calc.pack(side= tk.RIGHT,)

#create a label for the result or show when the input is wrong
labelr = tk.Label(text = "" , bd = 10)


#create the label for distance between the pots 
label_abstand = tk.Label(master = frame_abstand, text = "Distance between the pots", width = 20)
label_abstand.pack(side = tk.LEFT)

entry_abstand = tk.Entry(master = frame_abstand, width = 20)
entry_abstand.pack(side = tk.LEFT)

#pack all the frames, buttons and labels
frame_töpfe.pack( anchor="nw" )

frame_tiefe.pack(anchor="nw") 

frame_breite.pack(anchor="nw")

frame_abstand.pack(anchor="nw")

frame_knöpfe.pack(side = tk.BOTTOM, anchor="se")

labelr.pack(anchor = "nw")
    
# start the mainloop
window.mainloop()


>>>>>>> 6b558fa571cea4219737a2f491f4b55606cf1d58


