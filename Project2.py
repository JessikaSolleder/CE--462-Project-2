import tkinter as tk
from tkinter import messagebox
from tkinter import simpledialog
import numpy as np
import math
import matplotlib.pyplot as plt

# Access Euler's number & Pi  
euler_number = math.e
pi_value = math.pi

def main():

# # WALL GEOMETRY # #
    # Constants
    global GammaConcrete, StemTop, Depth, StemSlope
    GammaConcrete = 150 # pcf
    StemTop = 1 # ft
    Depth = 2 # ft
    StemSlope = 1/0.02 # unitless

    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a pop up window
    global Height, Gamma, Alpha, Phi
    Height = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")
    Gamma = simpledialog.askfloat("Input", "Enter the value of \u03B3 in lbs/cubic foot.")
    Alpha = simpledialog.askfloat("Input", "Enter the value of \u03B1 in degrees.")
    Phi = simpledialog.askfloat("Input", "Enter the value of \u03c6 ' in degrees.")

# # Calculating Dimensions of Rigid Retaining Wall Based on User Input of Height # #

    global StemHeight, SlabBottom, SlabThickness, Heel, Toe, SoilHeight
    StemHeight = 1 * Height
    SlabBottom = 0.5 * Height
    SlabThickness = 0.1 * Height
    Heel = 0.3 * Height
    Toe = 0.1 * Height
    SoilHeight = SlabThickness + StemHeight + Heel + np.tan(Alpha)


def rankine_analysis(SoilHeight, Gamma, Alpha, Phi):
    global RPv, RPh, RPa, RPa, RKp, RKa, RPp
    
    # Convert angles from degrees to radians
    
    global Phi_rad, Alpha_rad
    
    Phi_rad= np.radians(Phi)
    Alpha_rad= np.radians(Alpha)

    # Calculate the coefficient (Ka) of the active earth pressure using Rakine
    RKa = (1 - np.sin(Phi_rad)) / (1 + np.sin(Phi_rad))

    # Calculate the coefficient (Kp) of the passive earth pressure using Rakine
    RKp = (1 + np.sin(Phi_rad)) / (1 - np.sin(Phi_rad))

    # Calculate Rankine Lateral Earth Pressures for the active and passive conditions
    RPa = 0.5 * RKa * (Gamma * (SoilHeight ** 2))
    RPp = 0.5 * RKp * (Gamma * (SoilHeight ** 2))

    # Calculate Rankine horizontal & vertical components of the active earth pressure from above
    RPv = RPa * np.sin(Alpha_rad)
    RPh = RPa * np.cos(Alpha_rad)
    
#####################################################################################################################################
# ADD AN OUTPUT TABLE AND MAYBE AN IMAGE OR SOMETHING TO WRAP UP RANKINE
#####################################################################################################################################

def coloumb_pressure():
    global CPv, CPh, CPa, CPa, CKp, CKa, CPp

 # Calculate the coefficient (Ka) of the active earth pressure using Coulomb
    CKa = ((np.cos(Phi_rad)) ** 2) / (np.cos(Alpha_rad)*(1+((np.sin(Phi_rad + Alpha_rad) * (np.sin(Phi_rad - Alpha_rad)))/(np.cos(Alpha_rad)*(np.cos(-Alpha_rad))) ** 0.5) ** 2))

 # Calculate the coefficient (Kp) of the passive earth pressure using Coulomb
    CKp = (1 + np.sin(Phi_rad)) / (1 - np.sin(Phi_rad))
    
# Calculate Coulomb Lateral Earth Pressures for the active and passive conditions
    CPp = 0.5 * RKp * (Gamma * (SoilHeight ** 2))
    CPa = 0.5 * (CKa*Gamma*(SoilHeight ** 2))
    
# Calculate Coulomb horizontal & vertical components of the active earth pressure from above 

    CPv = CPa * np.sin(Alpha_rad)
    CPh = CPa * np.cos(Alpha_rad)
    
#####################################################################################################################################
# ADD AN OUTPUT TABLE AND MAYBE AN IMAGE OR SOMETHING TO WRAP UP COULOMB
#####################################################################################################################################

# # Overturning Calculations # #
# # See diagram included in GitHub Repository for diagram with sections labeled that will be referenced here # #
# # Calculating Moments Originating from Structure # #
# Take moment about the center point at the base where the stem meets with the slab, assume CCW = +
# Moment location (assuming the center of the bottom of the slab is (0,0)).
# (x,y) = (0.15 * H, 0.1 * H)

def overturning_moment (Height, Gamma, GammaConcrete, Heel, SoilHeight, SlabThickness, StemTop, Alpha_rad):
#####################################################################################################################################

# # Overturning Calculations # #
# # See diagram included in GitHub Repository for diagram with sections labeled that will be referenced here # #
# # Calculating Moments Originating from Structure # #
# Take moment about the center point at the base where the stem meets with the slab, assume CCW = +
# Moment location (assuming the center of the bottom of the slab is (0,0)).
# (x,y) = (0.15 * H, 0.1 * H)

    global WeightTotal, OverturnFSR, OverturnFSC, VerticalForceSum
    
    Area1 = ((Heel) * (Height - SlabThickness))
    Area2 = 0.5 * ((0.3 * Height) ** 2) * (np.tan(Alpha_rad))
    Area3 = StemTop * (Height - SlabThickness)
    Area4 = 0.5 * (Height * (0.02 * Height))
    Area5 = (0.1 * Height) * (0.5 * Height)

    Weight1 = Gamma * Area1
    Weight2 = Gamma * Area2
    Weight3 = GammaConcrete * Area3
    Weight4 = GammaConcrete * Area4 
    Weight5 = GammaConcrete * Area5
    WeightTotal = Weight1 + Weight2 + Weight3 + Weight4 + Weight5

    VerticalForce1 = Weight1
    VerticalForce2 = Weight2
    VerticalForce3 = Weight3
    VerticalForce4 = Weight4
    VerticalForce5 = Weight5
    VerticalForceSum = VerticalForce1 + VerticalForce2 + VerticalForce3 + VerticalForce4 + VerticalForce5


    MomentArm1 = (0.15 * Height)
    MomentArm2 = (2/3) * (0.3 * Height)
    MomentArm3 = ((0.1 * Height) - (0.02 * Height)) / 2
    MomentArm4 = 0.5 * ((0.1 * Height) - StemTop) ###### DOUBLE CHECK THIS ONE! #####
    MomentArm5 = (SlabBottom / 2) - (Toe + (0.05 * Height))

    Moment1 = VerticalForce1 * MomentArm1
    Moment2 = VerticalForce2 * MomentArm2
    Moment3 = VerticalForce3 * MomentArm3
    Moment4 = VerticalForce4 * MomentArm4
    Moment5 = VerticalForce5 * MomentArm5
    MomentSum =  - Moment1 - Moment2 - Moment3 + Moment4 - Moment5

    # Rankine: Resisting Moment
    MomentResistR = MomentSum + SlabBottom * RPv

    # Rankine: Overturning Moment
    MomentOverturnR = (SoilHeight / 3) * RPh

    # Coulomb: Resisting Moment
    MomentResistC = MomentSum + SlabBottom * CPv

    # Coulomb: Overturning Moment
    MomentOverturnC = (SoilHeight / 3) * CPh

    # Rankine: Overturning FS
    OverturnFSR = MomentResistR / MomentOverturnR

    # Coulomb: Overtruning FS
    OverturnFSC = MomentResistC / MomentOverturnC
    

def FS_Sliding():
#######################################################################
# FS Sliding
    global SlidingFSC, SlidingFSR
# Rankine: Sliding FS Check
    SigmaR = (3/4) * (Phi) # Interface friction between the concrete and the base soil
    SlidingFSR = (WeightTotal * np.tan(SigmaR)) / (RPh)

    # Coulomb: Sliding FS Check
    SigmaC = (3/4) * (Phi) # Interface friction between the concrete and the base soil
    SlidingFSC = (WeightTotal * np.tan(SigmaC)) / (CPh)
    
    
    
def qmin_qmax_qeq():
    ##############################################################################
    # Calculate qmin and qmax
    global qmax, qmin, qeq
    
    qmin = (VerticalForceSum / SlabBottom) * (1 - ((6 * np.exp)/SlabBottom)) #psf
    qmax = (VerticalForceSum / SlabBottom) * (1 + ((6 * np.exp)/SlabBottom)) #psf
    qeq =  (VerticalForceSum / Bprime) #psf

def Bearing_Capacity():
    ##############################################################################
    # Bearing Capacity, values will be the same whether user prefered to use Rankine 
    # or Coulomb above
    global Nq, Nc, Ngamma, q, Bprime, qult_terz, Bearing_Capacity_FS
    
    cprime = 0 # we are assuming there is no cohesion
    Nq = np.exp(np.pi * np.tan(Phi_rad)) * np.tan(np.radians(45) + Phi_rad / 2) ** 2
    Nc = (Nq - 1) * np.cot(Phi_rad)
    Ngamma = 2 * (Nq + 1) * np.tan(Phi_rad)
    
    q = Depth * Gamma #psf
    Bprime = SlabBottom - 2 * np.exp #ft
    qult_terz = (cprime * Nc) + (q * Nq) + (0.5 * Gamma * Bprime * Ngamma) #psf
    Bearing_Capacity_FS = (qult_terz / qmax) # elected to use the trapezoidal max stress (qult_terz / qmax)
    
def Schmertmann():
    
    global Es1, Es2, Es3, Es4, Hc1, Hc2, Hc3, Hc4
    
    p = qeq #psf
    Po = Depth * Gamma #psf
    Pop = (Depth + Be) * Gamma #psf
    Be = Bprime #ft
    DeltaP = p - Po #psf
    
    ##################### MESSAGE ABOUT SECTION AND WHAT WILL BE DONE, HERE ##########
    Es1 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 1 in tsf.") #tsf
    Es2 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 2 in tsf.") #tsf
    Es3 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 3 in tsf.") #tsf
    Es4 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 4 in tsf.") #tsf

    Hc1 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 1 in feet.") #ft
    Hc2 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 2 in feet.") #ft
    Hc3 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 3 in feet.") #ft
    Hc4 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 4 in feet.") #ft
    
    Iz1 = 0.5 + 0.1 * (DeltaP / Pop) ** 0.5
    


if __name__ == "__main__":
    main()