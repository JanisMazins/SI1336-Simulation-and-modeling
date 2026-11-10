#!/bin/python3

# Python simulation of a double pendulum with real time animation.
# BH, MP, AJ, TS 2020-10-27, latest version 2022-10-25.

from matplotlib import animation
from pylab import *
import numpy as np

"""
    This script simulates and animates a double pendulum.
    Classes are similar to the ones of pendolum_template.py. The main differences are:
    - coordinates are obtained from the total energy value E (look at what functions
        Oscillator.p2squaredFromH and Oscillator.__init__ do)
    - you are asked to implement the expression for the derivatives of the Hamiltonian 
        w.r.t. coordinates p1 and p2
    - you are asked to check when the conditions to produce the Poincare' map are
        satisfied and append the coordinates' values to some container to then plot
"""

# Global constants
G = 9.8  # gravitational acceleration

# Kinetic energy
def Ekin(osc):
    return 1 / (2.0 * osc.m * osc.L * osc.L) * (
            osc.p1 * osc.p1 + 2.0 * osc.p2 * osc.p2 - 2.0 * osc.p1 * osc.p2 * cos(osc.q1 - osc.q2)) / (
                   1 + (sin(osc.q1 - osc.q2)) ** 2)

# Potential energy
def Epot(osc):
    return osc.m * G * osc.L * (3 - 2 * np.cos(osc.q1) - np.cos(osc.q2))


# Class that holds the parameter and state of a double pendulum
class Oscillator:

    def p2squaredFromH(self):
        return (self.E - Epot(self)) * (1 + (sin(self.q1 - self.q2)) ** 2) * self.m * self.L * self.L

    # Initial condition is [q1, q2, p1, p2]; p2 is however re-obtained based on the value of E
    # therefore you can use any value for init_cond[3]
    def __init__(self, m=1, L=1, t0=0, E=15, init_cond=[0.0, 0.0, 0.0, -1.0]) :
        self.m = m      # mass of the pendulum bob
        self.L = L      # arm length
        self.t = t0     # the initial time
        self.E = E      # total conserved energy
        self.q1 = init_cond[0]
        self.q2 = init_cond[1]
        self.p1 = init_cond[2]
        self.p2 = -1.0
        while (self.p2 < 0):
            # Comment the two following lines in case you want to exactly prescribe values to q1 and q2
            # However, be sure that there exists a value of p2 compatible with the imposed total energy E!
            self.q1 = np.pi * (2 * np.random.random() - 1)
            self.q2 = np.pi * (2 * np.random.random() - 1)
            p2squared = self.p2squaredFromH()
            if (p2squared >= 0):
                self.p2 = np.sqrt(p2squared)
        self.q2_prev = self.q2
        print("Initialization:")
        print("E  = "+str(self.E))
        print("q1 = "+str(self.q1))
        print("q2 = "+str(self.q2))
        print("p1 = "+str(self.p1))
        print("p2 = "+str(self.p2))


# Class for storing observables for an oscillator
class Observables:

    def __init__(self):
        self.time = []          # list to store time
        self.q1list = []        # list to store q1
        self.q2list = []        # list to store q2
        self.p1list = []        # list to store p1
        self.p2list = []        # list to store p2
        self.epot = []          # list to store potential energy
        self.ekin = []          # list to store kinetic energy
        self.etot = []          # list to store total energy
        self.poincare_q1 = []   # list to store q1 for Poincare plot
        self.poincare_p1 = []   # list to store p1 for Poincare plot


# Derivate of H with respect to p1
def dHdp1(q1, q2, p1, p2, m, L):
    return (p1 - p2*np.cos(q1 - q2))/(L**2*m*((np.sin(q1 - q1))**2 + 1))


# Derivate of H with respect to p2
def dHdp2(q1, q2, p1, p2, m, L):
    return (2*p2 - p1*np.cos(q1 - q2))/(L**2*m*((np.sin(q1 - q2))**2 + 1))


# Derivate of H with respect to q1
def dHdq1(q1, q2, p1, p2, m, L):
    return 1 / (2.0 * m * L * L) * (
        -2 * (p1 * p1 + 2 * p2 * p2) * cos(q1 - q2) + p1 * p2 * (4 + 2 * (cos(q1 - q2)) ** 2)) * sin(
            q1 - q2) / (1 + (sin(q1 - q2)) ** 2) ** 2 + m * G * L * 2.0 * sin(q1)


# Derivate of H with respect to q2
def dHdq2(q1, q2, p1, p2, m, L):
    return 1 / (2.0 * m * L * L) * (
        2 * (p1 * p1 + 2 * p2 * p2) * cos(q1 - q2) - p1 * p2 * (4 + 2 * (cos(q1 - q2)) ** 2)) * sin(q1 - q2) / (
            1 + (sin(q1 - q2)) ** 2) ** 2 + m * G * L * sin(q2)


class RK4Integrator:

    def __init__(self, dt=0.01):
        self.dt = dt    # time step

    def integrate(self,
                  osc,
                  obs, 
                  ):

        """ Perform a single integration step """
        self.timestep(osc, obs)

        """ Append observables to their lists """
        obs.time.append(osc.t)
        obs.q1list.append(osc.q1)
        obs.q2list.append(osc.q2)
        obs.p1list.append(osc.p1)
        obs.p2list.append(osc.p2)
        obs.epot.append(Epot(osc))
        obs.ekin.append(Ekin(osc))
        obs.etot.append(Epot(osc) + Ekin(osc))
        if abs(osc.q2) < 10**-2 and osc.p2 - osc.p1 > 0:
            obs.poincare_q1.append(osc.q1)
            obs.poincare_p1.append(osc.p1)


    """
        Implementation of RK4 for a system of 4 variables
    """
    def timestep(self, osc, obs):
        dt = self.dt
        osc.t += dt

        # Initialization
        q1 = osc.q1
        q2 = osc.q2
        p1 = osc.p1
        p2 = osc.p2
        m = osc.m
        L = osc.L

        # RK4 coefficients (Butcher tableau)
        # Positions
        a = np.zeros( (2,4), dtype=float )
        # Momenta
        b = np.zeros( (2,4), dtype=float )

        # First sub-step:
        a[0,0] = dt * dHdp1(q1, q2, p1, p2, m, L)
        a[1,0] = dt * dHdp2(q1, q2, p1, p2, m, L)
        b[0,0] = - dt * dHdq1(q1, q2, p1, p2, m, L)
        b[1,0] = - dt * dHdq2(q1, q2, p1, p2, m, L)

        # Second sub-step:
        a[0,1] = dt * dHdp1(q1+0.5*a[0,0], q2+0.5*a[1,0], p1+0.5*b[0,0], p2+0.5*b[1,0], m, L) 
        a[1,1] = dt * dHdp2(q1+0.5*a[0,0], q2+0.5*a[1,0], p1+0.5*b[0,0], p2+0.5*b[1,0], m, L) 
        b[0,1] = - dt * dHdq1(q1+0.5*a[0,0], q2+0.5*a[1,0], p1+0.5*b[0,0], p2+0.5*b[1,0], m, L)
        b[1,1] = - dt * dHdq2(q1+0.5*a[0,0], q2+0.5*a[1,0], p1+0.5*b[0,0], p2+0.5*b[1,0], m, L) 

        # Third sub-step:
        a[0,2] = dt * dHdp1(q1+0.5*a[0,1], q2+0.5*a[1,1], p1+0.5*b[0,1], p2+0.5*b[1,1], m, L) 
        a[1,2] = dt * dHdp2(q1+0.5*a[0,1], q2+0.5*a[1,1], p1+0.5*b[0,1], p2+0.5*b[1,1], m, L) 
        b[0,2] = - dt * dHdq1(q1+0.5*a[0,1], q2+0.5*a[1,1], p1+0.5*b[0,1], p2+0.5*b[1,1], m, L)
        b[1,2] = - dt * dHdq2(q1+0.5*a[0,1], q2+0.5*a[1,1], p1+0.5*b[0,1], p2+0.5*b[1,1], m, L)

        # Fourth sub-step:
        a[0,3] = dt * dHdp1(q1+a[0,2], q2+a[1,2], p1+b[0,2], p2+b[1,2], m, L)
        a[1,3] = dt * dHdp2(q1+a[0,2], q2+a[1,2], p1+b[0,2], p2+b[1,2], m, L)
        b[0,3] = - dt * dHdq1(q1+a[0,2], q2+a[1,2], p1+b[0,2], p2+b[1,2], m, L)
        b[1,3] = - dt * dHdq2(q1+a[0,2], q2+a[1,2], p1+b[0,2], p2+b[1,2], m, L)

        osc.q1 += (1.0/6.0) * ( a[0,0]+2.0*a[0,1]+2.0*a[0,2]+a[0,3] )
        osc.q2 += (1.0/6.0) * ( a[1,0]+2.0*a[1,1]+2.0*a[1,2]+a[1,3] )
        osc.p1 += (1.0/6.0) * ( b[0,0]+2.0*b[0,1]+2.0*b[0,2]+b[0,3] )
        osc.p2 += (1.0/6.0) * ( b[1,0]+2.0*b[1,1]+2.0*b[1,2]+b[1,3] ) 


# Animation function which integrates a few steps and return a line for the pendulum
def animate(framenr, osc, obs, integrator, pendulum_lines, stepsperframe):
    for it in range(stepsperframe):
        integrator.integrate(osc, obs)

    x1 = np.sin(osc.q1)
    y1 = -np.cos(osc.q1)
    x2 = x1 + np.sin(osc.q2)
    y2 = y1 - np.cos(osc.q2)
    pendulum_lines.set_data([0, x1, x2], [0, y1, y2])
    return pendulum_lines,


class Simulation:

    def reset(self, osc=Oscillator()) :
        self.oscillator = osc
        self.obs = Observables()

    def __init__(self, osc=Oscillator()) :
        self.reset(osc)

    def run(self,
            integrator,
            tmax=30.,   # final time
            outfile='energy1.pdf'
            ):

        n = int(tmax / integrator.dt)

        for it in range(n):
            integrator.integrate(self.oscillator, self.obs)

    def run_animate(self,
            integrator,
            tmax=30.,           # final time
            stepsperframe=5,    # how many integration steps between visualising frames
            outfile='energy1.pdf'
            ):

        numframes = int(tmax / (stepsperframe * integrator.dt))

        plt.clf()

        # If you experience problems visualizing the animation try to comment/uncomment this line
        fig = plt.figure()

        ax = plt.subplot(xlim=(-2.2, 2.2), ylim=(-2.2, 2.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_lines, = ax.plot([], [], lw=5)

        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(fig, animate,  # init_func=init,
                                       fargs=[self.oscillator, self.obs, integrator, pendulum_lines, stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)
        
        # If you experience problems visualizing the animation try to comment/uncomment this line
        # plt.show()

        # If you experience problems visualizing the animation try to comment/uncomment this line 
        plt.waitforbuttonpress(10)

    # Plot coordinates and energies (to be called after running)
    def plot_observables(self, title="Double pendulum") :

        plt.figure()
        plt.title("p1-q1 phase plane plot")
        plt.xlabel('q1')
        plt.ylabel('p1')
        plt.plot(self.obs.q1list, self.obs.p1list)
        plt.tight_layout()  # adapt the plot area tot the text with larger fonts 

        plt.figure()
        plt.title("p2-q2 phase plane plot")
        plt.xlabel('q2')
        plt.ylabel('p2')
        plt.plot(self.obs.q2list, self.obs.p2list)
        plt.plot([0.0, 0.0], [min(self.obs.p2list), max(self.obs.p2list)], 'k--')
        plt.tight_layout()  # adapt the plot area tot the text with larger fonts 

        plt.figure()
        plt.title("Energy of the double pendulum")
        plt.xlabel('time')
        plt.ylabel('energy')
        plt.plot(self.obs.time, self.obs.epot, self.obs.time, self.obs.ekin, self.obs.time, self.obs.etot)
        plt.legend(('Epot', 'Ekin', 'Etot'))
        plt.tight_layout()  # adapt the plot area tot the text with larger fonts 

        plt.show()

def get_p2(p1, q1, q2, H, L=1, m=1):
    c1 = -p1*np.cos(q1 - q2)
    c2 = 1/2*(p1**2 - (2*m*L**2*(1 + (np.sin(q1 - q2))**2))*(H + m*G*L*(2*np.cos(q1) + np.cos(q2) - 3)))
    
    if (c1/2)**2 - c2 < 0:
        raise ValueError("No real p2 satisfy the given observables")
    else:
        a = -c1/2
        b = np.sqrt((c1/2)**2 - c2)
        print(f"p_1,2 = {a} ± {b}")
        return a + b, a - b

def plot_poincare(poincare_q1, poincare_p1):
    plt.figure()
    plt.title("p1-q1 Poincaré plot")
    plt.xlabel('q1')
    plt.ylabel('p1')
    plt.plot(poincare_q1, poincare_p1, 'ro', markersize=5)
    plt.tight_layout()  # adapt the plot area tot the text with larger fonts
    plt.show()

def exercise_15ab():
    t_max = 150
    dt = 0.003
    integrator = RK4Integrator()
    integrator.dt = dt

    #q1 = [0, 0.5, 1]*pi
    #q2 = [0, 0.5, 1]*pi
    #H = [1, 15, 40]
    q1 = 1.1058531550558046
    q2 = 3.061392570207309
    p1 = 0
    H = 12
    #p21, p22 = get_p2(p1, q1, q2, H)
    
    S = Simulation()
    S.reset(osc=Oscillator(m=1, L=1, t0=0, E=H))
    S.run_animate(integrator, tmax=t_max, stepsperframe=15)
    S.plot_observables()
    plot_poincare(S.obs.poincare_q1, S.obs.poincare_p1)

def exercise_15c():
    H = 15
    q1_12 = [0, 1.1]
    q2 = 0
    p1 = 0

    t_max = 2000
    dt = 0.003

    integrator = RK4Integrator()
    integrator.dt = dt

    poincare_q1_list = []
    poincare_p1_list = []
    for i in range(2):
        S = Simulation()
        q1 = q1_12[i]
        p2 = get_p2(p1, q1, q2, H)
        S.reset(osc=Oscillator(m=1, L=1, t0=0, E=H, init_cond=[q1, q2, p1, p2]))
        S.run(integrator, tmax=t_max)
        poincare_q1_list.extend(S.obs.poincare_q1)
        poincare_p1_list.extend(S.obs.poincare_p1)

    plot_poincare(poincare_q1_list, poincare_p1_list)

def exercise_15d():
    n = 10
    H = 25
    
    t_max = 500
    dt = 0.003

    integrator = RK4Integrator()
    integrator.dt = dt

    poincare_q1_list = []
    poincare_p1_list = []

    for i in range(n):
        print(f"Iteration: {i}")
        S = Simulation()
        S.reset(osc=Oscillator(m=1, L=1, t0=0, E=H))
        S.run(integrator, tmax=t_max)
        poincare_q1_list.extend(S.obs.poincare_q1)
        poincare_p1_list.extend(S.obs.poincare_p1)

    plot_poincare(poincare_q1_list, poincare_p1_list)





if __name__ == "__main__" :
    #get_p2(0, 0*np.pi, 0*np.pi, 1)
    exercise_15ab()
    #exercise_15c()
    #exercise_15d()
    # ...