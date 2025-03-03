# region imports
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
# endregion

# region functions
def ode_system(t, X, *params):
    '''
    The ode system is defined in terms of state variables.
    I have as unknowns:
    x: position of the piston (This is not strictly needed unless I want to know x(t))
    xdot: velocity of the piston
    p1: pressure on right of piston
    p2: pressure on left of the piston
    For initial conditions, we see: x=x0=0, xdot=0, p1=p1_0=p_a, p2=p2_0=p_a
    :param X: The list of state variables.
    :param t: The time for this instance of the function.
    :param params: the list of physical constants for the system.
    :return: The list of derivatives of the state variables.
    '''
    #unpack the parameters
    A, Cd, ps, pa, V, beta, rho, Kvalve, m, y=params

    #state variables
    #X[0]=x
    #X[1]=xdot
    #X[2]=p1
    #X[3]=p2

    #calculate derivitives
    #conveniently rename the state variables
    x = X[0]  # Position of the piston
    xdot = X[1]  # Velocity of the piston
    p1 = X[2]  # Pressure on the right side
    p2 = X[3]  # Pressure on the left side

    # Define the differential equations
    xddot = (p1 - p2) * A / m
    p1dot = (y * Kvalve * (ps - p1) - rho * A * xdot * beta / V)
    p2dot = (-y * Kvalve * (p2 - pa) - rho * A * xdot * beta / V)

    #return the list of derivatives of the state variables
    return [xdot, xddot, p1dot, p2dot]


def main():
    #After some trial and error, I found all the action seems to happen in the first 0.02 seconds
    t=np.linspace(0,0.02,200)
    #myargs=(A, Cd, Ps, Pa, V, beta, rho, Kvalve, m, y)
    myargs=(4.909E-4, 0.6, 1.4E7,1.0E5,1.473E-4,2.0E9,850.0,2.0E-5,30, 0.002)
    #because the solution calls for x, xdot, p1 and p2, I make these the state variables X[0], X[1], X[2], X[3]
    #ic=[x=0, xdot=0, p1=pa, p2=pa]
    pa = 1.0E5
    ic = [0, 0, myargs[3], myargs[3]]  # x = 0, xdot = 0, p1 = pa, p2 = pa
    #call odeint with ode_system as callback
    # Solve the system using solve_ivp
    t_span = (0, 0.02)
    t_eval = np.linspace(0, 0.02, 200)  # Time range
    sln=solve_ivp(ode_system, t_span, ic, args=myargs, t_eval=t_eval)

    #unpack result into meaningful names
    xvals=sln.y[0]
    xdot = sln.y[1]  # Velocity component from the ODE solution
    p1=sln.y[2]
    p2=sln.y[3]

    #plot the result
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.plot(sln.t, xvals, 'r-', label='$x$')
    plt.ylabel('$x$')
    plt.legend(loc='upper left')

    plt.twinx()
    plt.plot(sln.t, xdot, 'b-', label='$\dot{x}$ (Velocity)')
    plt.ylabel('$\dot{x}$ (m/s)')
    plt.legend(loc='lower right')

    # Plot p1 and p2
    plt.subplot(2, 1, 2)
    plt.plot(sln.t, p1, 'b-', label='$P_1$ (Pa)')
    plt.plot(sln.t, p2, 'r-', label='$P_2$ (Pa)')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure (Pa)')
    plt.legend(loc='lower right')

    plt.show()
# endregion

# region function calls
if __name__=="__main__":
    main()
# endregion