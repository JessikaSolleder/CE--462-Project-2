import tkinter as tk
from tkinter import messagebox
import numpy as np

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

# # Calculating Dimensions of Rigid Retaining Wall Based on User Input of Height # #
    global StemHeight, StemBottom, SlabThickness, Heel, Toe, SoilHeight
    StemHeight = 1 * Height
    StemBottom = 0.5 * Height
    SlabThickness = 0.1 * Height
    Heel = 0.1 * Height
    Toe = 0.3 * Height
    SoilHeight = SlabThickness + StemHeight + Heel + np.tan(Alpha)
 
def rankine_analysis(Height, Gamma, Alpha, Phi):
    # Recall the dimensions of the rigid retaining wall

    StemHeight = 1 * Height
    StemBottom = 0.5 * Height
    SlabThickness = 0.1 * Height
    Heel = 0.1 * Height
    Toe = 0.3 * Height
    SoilHeight = SlabThickness + StemHeight + Heel + np.tan(Alpha)
    
    # Convert angles from degrees to radians
    Phi_rad= np.radians(Phi)
    Alpha_rad= np.radians(Alpha)

    # Calculate the coefficient (Ka) of the active earth pressure using Rakine
    RKa = (1 - np.sin(Phi_rad)) / (1 + np.sin(Phi_rad))

    # Calculate Rankine passive earth pressure coefficient (Kp)
    RKp = (1 + np.sin(Phi_rad)) / (1 - np.sin(Phi_rad))

    # Calculate Rankine Lateral Earth Pressures for the active and passive conditions
    RPa = 0.5 * RKa * (Gamma * (SoilHeight ** 2))
    RPp = 0.5 * RKp * (Gamma * (SoilHeight ** 2))

    # Calculate Rankine horizontal & vertical components of the active earth pressure from above
    RPv = RPa * np.sin(Alpha_rad)
    RPh = RPa * np.cos(Alpha_rad)

def coloumb_pressure():

 # Calculate the coefficient (Ka) of the active earth pressure using Coulomb
 
    CKa = ((np.cos(Phi)) ** 2) / (np.cos(Alpha)*(1+((np.sin(Phi + Alpha) * (np.sin(Phi - Alpha)))/(np.cos(Alpha)*(np.cos(-Alpha))) ** 0.5) ** 2))

# Calculate Coulomb passive earth pressure coefficient (Kp)

# Calculate Coulomb Lateral Earth Pressures for the active and passive conditions

    CPa = 0.5 * (CKa*Gamma*(SoilHeight ** 2))
    CPv = CPa * np.sin(Alpha)
    
# Calculate Coulomb horizontal & vertical components of the active earth pressure from above 




if __name__ == "__main__":
    main()