# # WALL GEOMETRY # #
# GIVEN:

StemTop = 1
D = 2
StemSlope = 1/0.02

import tkinter as tk
from tkinter import simpledialog


def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a pop up window
    Height = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Height} , thank you!")

if __name__ == "__main__":
    main()
    


def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a popup window
    Gammab = simpledialog.askfloat("Input", "Enter the value of \u03B3b in lbs/cubic foot.")
    Gammat = simpledialog.askfloat("Input", "Enter the value of \u03B3t in lbs/ cubic foot.")
    Gammac = simpledialog.askfloat("Input", "Enter the value of \u03B3c in lbs/ cubic foot.")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Gammab, Gammat, Gammac} , for \u03B3b, \u03B3t and \u03B3c respectively, thank you!")

if __name__ == "__main__":
    main()
    
def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a popup window
    Phib = simpledialog.askfloat("Input", "Enter the value of \u03c6 'b in degrees")
    Phit = simpledialog.askfloat("Input", "Enter the value of \u03c6 't in degrees")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Phib, Phit} , for \u03c6 'b and \u03c6 't respectively, thank you!")

if __name__ == "__main__":
    main()

Stemheight = 1 * Height
Stembottom = 0.5 * Height
Slabthickness = 0.1 * Height
Heel = 0.1 * Height
Toe = 0.3 * Height


import numpy as np
import matplotlib.pyplot as plt

def rankine_pressure(depth, gamma, phi, delta):
    """
    Calculate lateral earth pressure using Rankine theory.

    Args:
    depth (numpy.ndarray): Array of depths of the soil surface from the ground surface (positive downwards) in meters.
    gamma (float): Unit weight of soil in kN/m^3.
    phi (float): Angle of internal friction of soil in degrees.
    delta (float): Angle of wall inclination in degrees.

    Returns:
    numpy.ndarray: Array of lateral earth pressure in kN/m^2 corresponding to each depth.
    """
    # Convert angles from degrees to radians
    phi = np.radians(phi)
    delta = np.radians(delta)

    # Calculate lateral earth pressure
    K = 1 - np.sin(phi)
    pressure = gamma * depth * K * np.cos(delta)

    return pressure

def plot_pressure(depth, pressure):
    
    # # Plot lateral earth pressure distribution.

    # Args:
    # depth (numpy.ndarray): Array of depths of the soil surface from the ground surface (positive downwards) in meters.
    # pressure (numpy.ndarray): Array of lateral earth pressure in kN/m^2 corresponding to each depth.
    
    plt.figure(figsize=(8, 6))
    plt.plot(pressure, depth, color='blue', linewidth=2)
    plt.xlabel('Lateral Earth Pressure (kN/m²)')
    plt.ylabel('Depth (m)')
    plt.title('Lateral Earth Pressure Distribution')
    plt.grid(True)
    plt.gca().invert_yaxis()  # Invert y-axis to show depth increasing downwards
    plt.show()

def main():
    # Define parameters
    depth = np.linspace(0, 10, 100)  # Depths from 0 to 10 meters
    gamma = 18.5  # Unit weight of soil (kN/m^3)
    phi = 30  # Angle of internal friction of soil (degrees)
    
import math

slope = 1 / 0.02
delta = math.atan(slope)
delta_degrees = math.degrees(delta)
print("Angle of inclination (delta) in radians:", delta)
print("Angle of inclination (delta) in degrees:", delta_degrees)
# Angle of wall inclination (degrees)

# Calculate lateral earth pressure using Rankine theory
lateral_pressure = rankine_pressure(depth, gamma, phi, delta)

    # Plot lateral earth pressure distribution
plot_pressure(depth, lateral_pressure)

if __name__ == "__main__":
    main()
