import tkinter as tk
from tkinter import simpledialog
import numpy as np
import math
import matplotlib.pyplot as plt

# Constants
GAMMA_CONCRETE = 150  # pcf

def get_user_input():
    """
    Prompt the user to input dimensions using a pop-up window and store global variables.
    """
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    height = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")
    gamma = simpledialog.askfloat("Input", "Enter the value of \u03B3 in lbs/cubic foot.")
    alpha = simpledialog.askfloat("Input", "Enter the value of \u03B1 in degrees.")
    phi = simpledialog.askfloat("Input", "Enter the value of \u03c6 ' in degrees.")

    return height, gamma, alpha, phi

def calculate_wall_dimensions(height):
    """
    Calculate and return dimensions of the rigid retaining wall based on user input height.
    """
    stem_height = 1 * height
    slab_bottom = 0.5 * height
    slab_thickness = 0.1 * height
    heel = 0.3 * height
    toe = 0.1 * height
    soil_height = slab_thickness + stem_height + heel + np.tan(np.radians(alpha))
    
    return stem_height, slab_bottom, slab_thickness, heel, toe, soil_height

def rankine_analysis(soil_height, gamma, alpha, phi):
    """
    Perform Rankine analysis and plot lateral earth pressure distribution.
    """
    phi_rad = np.radians(phi)
    alpha_rad = np.radians(alpha)
    rka = (1 - np.sin(phi_rad)) / (1 + np.sin(phi_rad))
    rkp = (1 + np.sin(phi_rad)) / (1 - np.sin(phi_rad))

    rpa = 0.5 * rka * (gamma * (soil_height ** 2))
    rpp = 0.5 * rkp * (gamma * (soil_height ** 2))
    rpv = rpa * np.sin(alpha_rad)
    rph = rpa * np.cos(alpha_rad)

    depth = np.linspace(0, soil_height, 100)  # Generate 100 points from top to bottom of the wall
    pressure_active = rka * gamma * depth
    pressure_passive = rkp * gamma * depth

    plot_pressure_distribution(depth, pressure_active, pressure_passive)

def plot_pressure_distribution(depth, pressure_active, pressure_passive):
    """
    Plot the lateral earth pressure distribution.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(pressure_active, depth, label='Active Earth Pressure', color='blue')
    plt.plot(pressure_passive, depth, label='Passive Earth Pressure', color='red')
    plt.gca().invert_yaxis()
    plt.xlabel('Pressure (pcf)')
    plt.ylabel('Depth (ft)')
    plt.title('Lateral Earth Pressure Distribution')
    plt.legend()
    plt.grid(True)
    plt.show()

# Add your additional functions here (Coulomb analysis, Overturning calculations, etc.)

if __name__ == "__main__":
    height, gamma, alpha, phi = get_user_input()
    dimensions = calculate_wall_dimensions(height)
    rankine_analysis(*dimensions[4:], gamma, alpha, phi)
    # Call other analysis functions as needed, passing appropriate arguments
