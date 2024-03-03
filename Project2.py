# # WALL GEOMETRY # #
# Constant Dimensions (in Feet):

StemTop = 1
D = 2
StemSlope = 1/0.02

import tkinter as tk
from tkinter import simpledialog


def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a pop up window
    Heightinput = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Heightinput}!")

if __name__ == "__main__":
    main()

def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a popup window
    Gamma = simpledialog.askfloat("Input", "Enter the value of \u03B3b in lbs/cubic foot.")
    Alpha = simpledialog.askfloat("Input", "Enter the value of \u03B1 in degrees." )
    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Gamma } for \u03B3 !")

if __name__ == "__main__":
    main()
    
def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a popup window
    Phi = simpledialog.askfloat("Input", "Enter the value of \u03c6 ' in degrees.")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Phi} for \u03c6 !")

if __name__ == "__main__":
    main()
    
# # Calculating Dimensions of Rigid Retaining Wall Based on User Input of Height # #

Stemheight = 1 * Height
Stembottom = 0.5 * Height
Slabthickness = 0.1 * Height
Heel = 0.1 * Height
Toe = 0.3 * Height
Soilheight = Slabthickness + Stemheight + Heel + tan (Alpha)

import numpy as np
import matplotlib.pyplot as plt

def rankine_coefficients:

    # Convert angles from degrees to radians
    Phi = np.radians(Phi)
    Alpha = np.radians(Alpha)

    # Calculate active earth pressure using Rakine
    Rka = (1 - np.sin(Phi)) / (1 + np.sin(Phi)

    # Calculate Rankine active and passive coefficients
    Rkp = (1 + np.sin(Phi)) / (1 - np.sin(Phi))

    # Calculate Rankine Lateral Earth Pressure
    Rpa = 0.5 * Ka * ( Gamma * (Height ** 2))
    Rpp = 0.5 * Kp * ( Gamma * (Height ** 2))

    # Calculate Rankine horizontal & vertical components of the active earth pressure from above
    Rpv = Rpa * np.sin(Alpha)
    Rph = Rpa * np.cos(Alpha)

def coloumb_pressure:

 # Calculate active earth pressure using Coulomb
    Rka = (1 - np.sin(Phi)) / (1 + np.sin(Phi)

    # Calculate Rankine active and passive coefficients
    Rkp = (1 + np.sin(Phi)) / (1 - np.sin(Phi))

    # Calculate Rankine Lateral Earth Pressure
    Rpa = 0.5 * Ka * ( Gamma * (Height ** 2))
    Rpp = 0.5 * Kp * ( Gamma * (Height ** 2))

    # Calculate Rankine Horizontal & Vertical Components of the Active Earth Pressure from Above
    Rpv = Rpa * np.sin(Alpha)
    Rph = Rpa * np.cos(Alpha)

