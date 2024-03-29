import tkinter as tk
from tkinter import messagebox
from tkinter import simpledialog
import numpy as np
import math
import matplotlib.pyplot as plt

# Constants / numbers
euler_number = math.e
pi_value = math.pi
gamma_concrete = 150 # pcf
stem_top = 1 # ft
slab_depth = 2 # ft
stem_slope = 1/0.02 # unitless

def main_program():
    
    user_input()
    rankine_analysis(soil_height, gamma, alpha, phi)
    coloumb_analysis()
    # Add other function calls and print statements as necessary
    overturning_moment(height, gamma, gamma_concrete, heel, soil_height, slab_thickness, stem_top, alpha_rad)
    FS_Sliding()
    qmin_qmax_qeq()
    Bearing_Capacity()
    Schmertmann_Method()

def user_input():
 
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    # Prompt the user to input dimensions using a pop up window
    
    global height, gamma, alpha, phi, alpha_rad
    height = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")
    gamma = simpledialog.askfloat("Input", "Enter the value of \u03B3 in lbs/cubic foot.")
    alpha = simpledialog.askfloat("Input", "Enter the value of \u03B1 in degrees.")
    phi = simpledialog.askfloat("Input", "Enter the value of \u03c6 ' in degrees.")

    # Calculating Dimensions of Rigid Retaining Wall Based on User Input of Height 
    
    global stem_height, slab_bottom, slab_thickness, heel, toe, soil_height
 
    alpha_rad= np.radians(alpha)
    stem_height = 1 * height
    slab_bottom = 0.5 * height
    slab_thickness = 0.1 * height
    heel = 0.3 * height
    toe = 0.1 * height
    soil_height = slab_thickness + stem_height + heel + np.tan(alpha_rad)


def rankine_analysis(soil_height, gamma, alpha, phi):
    global RPv, RPh, RPa, RPa, RKp, RKa, RPp, phi_rad, alpha_rad
    
    phi_rad= np.radians(phi)
    alpha_rad= np.radians(alpha)

    # Calculate the coefficient (Ka) of the active earth pressure using Rakine
    RKa = (1 - np.sin(phi_rad)) / (1 + np.sin(phi_rad))

    # Calculate the coefficient (Kp) of the passive earth pressure using Rakine
    RKp = (1 + np.sin(phi_rad)) / (1 - np.sin(phi_rad))

    # Calculate Rankine Lateral Earth Pressures for the active and passive conditions
    RPa = 0.5 * RKa * (gamma * (soil_height ** 2))
    RPp = 0.5 * RKp * (gamma * (soil_height ** 2))

    # Calculate Rankine horizontal & vertical components of the active earth pressure from above
    RPv = RPa * np.sin(alpha_rad)
    RPh = RPa * np.cos(alpha_rad)
    
    depth = np.linspace(0, soil_height, 100)  # Generate 100 points from top to bottom of the wall
    pressure_active = RKa * gamma * depth  # Active earth pressure distribution
    pressure_passive = RKp * gamma * depth  # Passive earth pressure distribution

    plt.figure(figsize=(10, 6))
    plt.plot(pressure_active, depth, label='Active Earth Pressure', color='blue')
    plt.plot(pressure_passive, depth, label='Passive Earth Pressure', color='red')
    plt.gca().invert_yaxis()  # Invert y-axis to show depth increasing downwards
    plt.xlabel('Pressure (pcf)')
    plt.ylabel('Depth (ft)')
    plt.title('Lateral Earth Pressure Distribution- Rankine')
    plt.legend()
    plt.grid(True)
    plt.show()
    
#####################################################################################################################################
# ADD AN OUTPUT TABLE AND MAYBE AN IMAGE OR SOMETHING TO WRAP UP RANKINE
#####################################################################################################################################

def coloumb_analysis():
    global CPv, CPh, CPa, CPa, CKp, CKa, CPp

 # Calculate the coefficient (Ka) of the active earth pressure using Coulomb
    CKa = ((np.cos(phi_rad)) ** 2) / (np.cos(alpha_rad)*(1+((np.sin(phi_rad + alpha_rad) * (np.sin(phi_rad - alpha_rad)))/(np.cos(alpha_rad)*(np.cos( - alpha_rad))) ** 0.5) ** 2))

 # Calculate the coefficient (Kp) of the passive earth pressure using Coulomb
    CKp = (1 + np.sin(phi_rad)) / (1 - np.sin(phi_rad))
    
# Calculate Coulomb Lateral Earth Pressures for the active and passive conditions
    CPp = 0.5 * RKp * (gamma * (soil_height ** 2))
    CPa = 0.5 * (CKa * gamma * (soil_height ** 2))
    
# Calculate Coulomb horizontal & vertical components of the active earth pressure from above 

    CPv = CPa * np.sin(alpha_rad)
    CPh = CPa * np.cos(alpha_rad)

#####################################################################################################################################
# ADD AN OUTPUT TABLE AND MAYBE AN IMAGE OR SOMETHING TO WRAP UP COULOMB
#####################################################################################################################################

# # Overturning Calculations # #
# # See diagram included in GitHub Repository for diagram with sections labeled that will be referenced here # #
# # Calculating Moments Originating from Structure # #
# Take moment about the center point at the base where the stem meets with the slab, assume CCW = +
# Moment location (assuming the center of the bottom of the slab is (0,0)).
# (x,y) = (0.15 * H, 0.1 * H)

def overturning_moment (height, gamma, gamma_concrete, heel, soil_height, slab_thickness, stem_top, alpha_rad):
#####################################################################################################################################

    global weight_total, overturn_FSR, overturn_FSC, vertical_force_sum, moment_resist_R, moment_overturn_R, moment_resist_C, moment_overturn_C
    
    area_1 = ((heel) * (height - slab_thickness))
    area_2 = 0.5 * ((0.3 * height) ** 2) * (np.tan(alpha_rad))
    area_3 = stem_top * (height - slab_thickness)
    area_4 = 0.5 * (height * (0.02 * height))
    area_5 = (0.1 * height) * (0.5 * height)

    weight_1 = gamma * area_1
    weight_2 = gamma * area_2
    weight_3 = gamma_concrete * area_3
    weight_4 = gamma_concrete * area_4 
    weight_5 = gamma_concrete * area_5
    weight_total = weight_1 + weight_2 + weight_3 + weight_4 + weight_5

    vertical_force_1 = weight_1
    vertical_force_2 = weight_2
    vertical_force_3 = weight_3
    vertical_force_4 = weight_4
    vertical_force_5 = weight_5
    vertical_force_sum = vertical_force_1 + vertical_force_2 + vertical_force_3 + vertical_force_4 + vertical_force_5


    moment_arm_1 = (0.15 * height)
    moment_arm_2 = (2/3) * (0.3 * height)
    moment_arm_3 = ((0.1 * height) - (0.02 * height)) / 2
    moment_arm_4 = 0.5 * ((0.1 * height) - stem_top) ###### DOUBLE CHECK THIS ONE! #####
    moment_arm_5 = (slab_bottom / 2) - (toe + (0.05 * height))

    moment_1 = vertical_force_1 * moment_arm_1
    moment_2 = vertical_force_2 * moment_arm_2
    moment_3 = vertical_force_3 * moment_arm_3
    moment_4 = vertical_force_4 * moment_arm_4
    moment_5 = vertical_force_5 * moment_arm_5
    moment_sum =  - moment_1 - moment_2 - moment_3 + moment_4 - moment_5

    # Rankine: Resisting Moment
    moment_resist_R = moment_sum + slab_bottom * RPv

    # Rankine: Overturning Moment
    moment_overturn_R = (soil_height / 3) * RPh

    # Coulomb: Resisting Moment
    moment_resist_C = moment_sum + slab_bottom * CPv

    # Coulomb: Overturning Moment
    moment_overturn_C = (soil_height / 3) * CPh

    # Rankine: Overturning FS
    overturn_FSR = moment_resist_R / moment_overturn_R

    # Coulomb: Overtruning FS
    overturn_FSC = moment_resist_C / moment_overturn_C
    

def FS_Sliding():
#######################################################################
# FS Sliding
    global sliding_FSC, sliding_FSR
    
# Rankine: Sliding FS Check
    sigma_R = (3/4) * (phi) # Interface friction between the concrete and the base soil
    sliding_FSR = (weight_total * np.tan(sigma_R)) / (RPh)

    # Coulomb: Sliding FS Check
    sigma_C = (3/4) * (phi) # Interface friction between the concrete and the base soil
    sliding_FSC = (weight_total * np.tan(sigma_C)) / (CPh)
    
    
    
def qmin_qmax_qeq():
    ##############################################################################
    # Calculate qmin and qmax
    global q_max, q_min, q_eq
    
    q_min = (vertical_force_sum / slab_bottom) * (1 - ((6 * math.e)/slab_bottom))  # psf
    q_max = (vertical_force_sum / slab_bottom) * (1 + ((6 * math.e)/slab_bottom))  # psf
    q_eq =  (vertical_force_sum / B_prime) #psf

def Bearing_Capacity():
    ##############################################################################
    # Bearing Capacity, values will be the same whether user prefered to use Rankine 
    # or Coulomb above
    global Nq, Nc, N_gamma, q, B_prime, q_ult_terz, bearing_capacity_FS, c_prime
    
    c_prime = 0 # we are assuming there is no cohesion
    Nq = np.exp(np.pi * np.tan(phi_rad)) * np.tan(np.radians(45) + phi_rad / 2) ** 2
    Nc = (Nq - 1) * (np.cos(phi_rad)/np.sin(phi_rad))
    N_gamma = 2 * (Nq + 1) * np.tan(phi_rad)
    
    q = slab_depth * gamma #psf
    B_prime = slab_bottom - 2 * math.e #ft
    q_ult_terz = (c_prime * Nc) + (q * Nq) + (0.5 * gamma * B_prime * N_gamma) #psf
    bearing_capacity_FS = (q_ult_terz / (vertical_force_sum / slab_bottom) * (1 + ((6 * math.e)/slab_bottom)) ) # elected to use the trapezoidal max stress (qult_terz / qmax)
    
    
def Schmertmann_Method():
    
    global Es_1, Es_2, Es_3, Es_4, Hc_1, Hc_2, Hc_3, Hc_4, Iz, t, c_1, c_2, delta_Hi_1, delta_Hi_2, delta_Hi_3, delta_Hi_4, delta_Hi_sum,S_i
    
    B_e = B_prime #ft
    p = q_eq #psf
    p_o = slab_depth * gamma #psf
    p_op = (slab_depth + B_e) * gamma #psf

    delta_p = p - p_o #psf
    
    ##################### MESSAGE ABOUT SECTION AND WHAT WILL BE DONE, HERE ##########
    Es_1 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 1 in tsf.") #tsf
    Es_2 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 2 in tsf.") #tsf
    Es_3 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 3 in tsf.") #tsf
    Es_4 = simpledialog.askfloat("Input", "Enter the modulus of elasticity for soil layer 4 in tsf.") #tsf

    Hc_1 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 1 in feet.") #ft
    Hc_2 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 2 in feet.") #ft
    Hc_3 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 3 in feet.") #ft
    Hc_4 = simpledialog.askfloat("Input", "Enter the thickness of soil layer 4 in feet.") #ft
    
    Iz = 0.5 + 0.1 * (delta_p / p_op) ** 0.5
    
    t = simpledialog.askfloat("Input", "Enter the number of years you want the settlement to be evaluated for.")
    
    delta_Hi_1 = Hc_1 * (Iz / 1 * Es_1)
    delta_Hi_2 = Hc_2 * (Iz / 1 * Es_2)
    delta_Hi_3 = Hc_3 * (Iz / 1 * Es_3)
    delta_Hi_4 = Hc_4 * (Iz / 1 * Es_4)
    delta_Hi_sum = delta_Hi_1 + delta_Hi_2 + delta_Hi_3 + delta_Hi_4
    
    c_1 = 1 - 0.5 * (p_o / delta_p) # >=0.5
    c_2 = 1 + 0.2 * np.log10(t/0.1)
    
    S_i = (c_1 * c_2 * delta_p) * delta_Hi_sum #feet
    
if __name__ == "__main__":
    user_input()
    # Temporarily remove to see the plot window
    # root.withdraw()  

    rankine_analysis(soil_height, gamma, alpha, phi)
    
    coloumb_analysis()
    print("Rankine: Vertical Component of Active Earth Pressure: ", RPv)
    print("Rankine: Horizontal Component of Active Earth Pressure: ", RPh)
    print("Coulomb: Vertical Component of Active Earth Pressure:", CPv)
    print("Coulomb: Horizontal Component of Active Earth Pressure:", CPh)
    
    overturning_moment (height, gamma, gamma_concrete, heel, soil_height, slab_thickness, stem_top, alpha_rad)
    
    print()
    print("Rankine: Total Resisting Moment: ", moment_resist_R)
    print("Rankine: Total Overturning Moment: ", moment_overturn_R)
    print("Rankine: Overturning Moment FS: ", overturn_FSR)
    print("Coulomb: Total Resisting Moment: ", moment_resist_C)
    print("Coulomb: Total Overturning Moment: ", moment_overturn_C)
    print("Coulomb: Overturning Moment FS: ", overturn_FSC)
    
    Bearing_Capacity()
    
    print()
    print("This program assumes that the soil is cohesionless. Using Terzaghi's method, the following bearing capacity information was calculated:")
    print("Value of Nq: ", Nq)
    print("Value of Nc: ", Nc)
    print("Value of q: ", q)
    print("Value of B': ", B_prime)
    print("Value of q_ultimate: ", q_ult_terz)
    print("Bearing Capacity FS: ", bearing_capacity_FS)
    
    FS_Sliding()
    
    print()
    print("Rankine: Sliding FS: ", sliding_FSR)
    print("Coulomb: Sliding FS: ", sliding_FSC)
    qmin_qmax_qeq()
    
    Schmertmann_Method()
    
    print()
    print("Schmertmann's Method: C1: ", c_1)
    print("Schmertmann's Method: C2: ", c_2)
    print("Schmertmann's Method: Immediate Settlement in Feet ", S_i)

 
continue_program = True
while continue_program:
        main_program()
        root = tk.Tk()
        root.withdraw()  # Hide the root window
        # Ask user if they want to rerun the program
        user_choice = messagebox.askyesno("Continue", "Would you like to rerun the program with new inputs?")
        if not user_choice:
            continue_program = False
        root.destroy()
        