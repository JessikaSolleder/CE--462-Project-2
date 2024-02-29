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

    # Prompt the user to input dimensions using a popup window
    Height = simpledialog.askfloat("Input", "Enter the height of the stem (H) in feet.")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Height} , thank you!")

if __name__ == "__main__":
    main()
    
def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a popup window
    Gammab = simpledialog.askfloat("Input", "Enter the value of Gamma b")
    Gammat = simpledialog.askfloat("Input", "Enter the value of Gamma t")
    Gammac = simpledialog.askfloat("Input", "Enter the value of Gamma c")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Gammab, Gammat, Gammac} , for Gammab, Gammat and Gammac respectively, thank you!")

if __name__ == "__main__":
    main()
    
def main():
    root = tk.Tk()
    root.withdraw()  # Hide the root window

    # Prompt the user to input dimensions using a popup window
    Phib = simpledialog.askfloat("Input", "Enter the value of Phi sub b")
    Phit = simpledialog.askfloat("Input", "Enter the value of Phi sub t")

    # Display the result in a message box
    tk.messagebox.showinfo("Confirmed", f"You input: {Phib, Phit} , for Phi sub b and Phi sub t respectively, thank you!")

if __name__ == "__main__":
    main()

## Did it work?##
# hOW ABIUT NOW?
# whagt about now ven?