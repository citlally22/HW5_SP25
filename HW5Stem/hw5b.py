# region imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import hw5a as pta
import random as rnd
from matplotlib import pyplot as plt
# endregion

# Constants
GALLON_TO_CUBIC_FEET = 0.133681  # 1 gallon = 0.133681 cubic feet
INCH_TO_FEET = 1/12  # 1 inch = 1/12 feet
G = 32.174  # ft/s^2, acceleration due to gravity

# Function to compute Reynolds number
def reynolds_number(diameter, velocity, kinematic_viscosity=1.1e-5):
    return (velocity * diameter) / kinematic_viscosity

# Function to compute head loss per foot
def head_loss(f, L=1, d=1, V=1):
    return (f * (L / d) * (V ** 2) / (2 * G))

# Function to calculate friction factor with transition region randomness
def ffPoint(Re, rr):
    """
    This function takes Re and rr as parameters and outputs a friction factor according to the following:
    1.  if Re>4000 use Colebrook Equation
    2.  if Re<2000 use f=64/Re
    3.  else calculate a probabilistic friction factor where the distribution has a mean midway between the prediction
        of the f=64/Re and Colebrook Equations and a standard deviation of 20% of this mean
    :param Re:  the Reynolds number
    :param rr:  the relative roughness
    :return:  the friction factor
    """
    if Re>=4000:
        return pta.ff(Re, rr,CBEQN=True)
    if Re<=2000:
        return pta.ff(Re, rr)

    #Compute friction factors for transition region
    CBff= pta.ff(Re, rr, CBEQN=True)  # Colebrook Equation prediction
    Lamff= pta.ff(Re, 0)  # Laminar Equation prediction (using smooth pipe assumption)
    mean=(CBff+Lamff)/2 # Midway mean value
    sig=0.2*mean # 20% standard deviation
    return rnd.normalvariate(mean, sig)  # Randomized friction factor for transition region

# Function to plot the computed point on the Moody diagram
def PlotPoint(Re, f, transition=False):
    pta.plotMoody(plotPoint=True, pt=(Re, f)) #call plotMoody without 'marker'
    print(f"DEBUG: Re={Re}, f={f}, transition={transition}")  # Debugging output

    marker = '^' if transition else 'o'  # Triangle for transition, circle otherwise
    print(f"DEBUG: Using marker '{marker}' for plotting")  # Confirm marker selection

    plt.plot(Re, f, marker, markersize=8, markeredgecolor='red', markerfacecolor='none')
    plt.show(block=False)
    plt.pause(0.1)


# Main function to get user input and compute friction factor
def main():
    while True:
        try:
            print("\n--- Enter Pipe Parameters ---")  # Clear separation for each input set
            d_inch = float(input("Enter pipe diameter (in inches): "))
            epsilon_mic = float(input("Enter pipe roughness (micro-inches): "))
            Q_gpm = float(input("Enter flow rate (gallons per minute): "))

            # Convert units
            d = d_inch * INCH_TO_FEET  # Convert inches to feet
            epsilon = epsilon_mic * 1e-6 * INCH_TO_FEET  # Convert micro-inches to feet
            Q = Q_gpm * GALLON_TO_CUBIC_FEET  # Convert GPM to cubic feet per second
            A = np.pi * (d / 2) ** 2  # Cross-sectional area
            V = Q / A  # Velocity
            Re = reynolds_number(d, V)  # Compute Reynolds number
            rr = epsilon / d  # Relative roughness
            f = ffPoint(Re, rr)  # Compute friction factor
            hf_per_L = head_loss(f, d=d, V=V)  # Compute head loss per foot

            print(f"\nReynolds number calculated: {Re}")
            print(f"Computed Friction Factor: {f:.5f}")
            print(f"Head Loss per foot: {hf_per_L:.5f} ft/ft")

            # Determine transition condition
            if 2000 < Re < 4000:
                transition = True
            else:
                transition = False

            # Plot the computed point
            PlotPoint(Re, f, transition)

            # Ask user if they want to input another set
            cont = input("\nDo you want to enter another set? (y/n): ").strip().lower()
            print(f"You entered: {cont}")  # Debugging line to check input

            if cont not in ['y', 'yes']:
                print("Exiting program.")
                break  # Exit the loop if the user enters anything other than 'y' or 'yes'

        except ValueError:
            print("\nInvalid input. Please enter numerical values.")
            continue  # Ensure the loop continues after an invalid input

if __name__ == "__main__":
    main()
