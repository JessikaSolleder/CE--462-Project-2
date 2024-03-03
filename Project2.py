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
    global Height, Gamma, Alpha, Phi
    Height = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")
    Gamma = simpledialog.askfloat("Input", "Enter the value of \u03B3 in lbs/cubic foot.")
    Alpha = simpledialog.askfloat("Input", "Enter the value of \u03B1 in degrees.")
    Phi = simpledialog.askfloat("Input", "Enter the value of \u03c6 ' in degrees.")


    Stemheight = 1 * Height
    Stembottom = 0.5 * Height
    Slabthickness = 0.1 * Height
    Heel = 0.1 * Height
    Toe = 0.3 * Height
    Height = Slabthickness + Stemheight + Heel + np.tan * (Alpha)
 
    
# # Calculating Dimensions of Rigid Retaining Wall Based on User Input of Height # #

import numpy as np
import matplotlib.pyplot as plt

def rankine_coefficients():

    # Convert angles from degrees to radians
    Phi = np.radians(Phi)
    Alpha = np.radians(Alpha)

    Stemheight = 1 * Height
    Stembottom = 0.5 * Height
    Slabthickness = 0.1 * Height
    Heel = 0.1 * Height
    Toe = 0.3 * Height
    Height = Slabthickness + Stemheight + Heel + np.tan * (Alpha)
    

    # Calculate the coefficient (Ka) of the active earth pressure using Rakine
    Rka = (1 - np.sin(Phi)) / (1 + np.sin(Phi))

    # Calculate Rankine passive earth pressure coefficient (Kp)
    Rkp = (1 + np.sin(Phi)) / (1 - np.sin(Phi))

    # Calculate Rankine Lateral Earth Pressures for the active and passive conditions
    Rpa = 0.5 * Rka * ( Gamma * (Height ** 2))
    Rpp = 0.5 * Rkp * ( Gamma * (Height ** 2))

    # Calculate Rankine horizontal & vertical components of the active earth pressure from above
    Rpv = Rpa * np.sin(Alpha)
    Rph = Rpa * np.cos(Alpha)
    
    
def coloumb_pressure():

 # Calculate the coefficient (Ka) of the active earth pressure using Coulomb
 
    Cka = ((np.cos(Phi)) ** 2) / (np.cos(Alpha)*(1+((np.sin(Phi + Alpha) * (np.sin(Phi - Alpha)))/(np.cos(Alpha)*(np.cos(-Alpha))) ** 0.5) ** 2))

# Calculate Coulomb passive earth pressure coefficient (Kp)

# Calculate Coulomb Lateral Earth Pressures for the active and passive conditions

# Calculate Coulomb horizontal & vertical components of the active earth pressure from above 


if __name__ == "__main__":
    main()