#!/bin/python3

# Python simulation of a simple planar pendulum with real time animation
# BH, OF, MP, AJ, TS 2020-10-20, latest version 2022-10-25.

from matplotlib import animation
from pylab import *

# Global constants
G = 9.8  # gravitational acceleration

class Oscillator:

    """ Class for a general, simple oscillator """

    def __init__(self, m=1, c=9, t0=0, theta0=np.pi*0.5, dtheta0=0, gamma=0):
        self.m = m              # mass of the pendulum bob
        self.c = c              # c = g/L
        self.L = G / c          # string length
        self.t = t0             # the time
        self.theta = theta0     # the position/angle
        self.dtheta = dtheta0   # the velocity
        self.gamma = gamma      # damping

class Observables:

    """ Class for storing observables for an oscillator """

    def __init__(self):
        self.time = []          # list to store time
        self.pos = []           # list to store positions
        self.vel = []           # list to store velocities
        self.energy = []        # list to store energy


class BaseSystem:
    
    def force(self, osc):

        """ Virtual method: implemented by the child classes  """

        pass


class Harmonic(BaseSystem):
    def force(self, osc):
        return - osc.m * ( osc.c*osc.theta + osc.gamma*osc.dtheta )


class Pendulum(BaseSystem):
    def force(self, osc):
        return - osc.m * ( osc.c*np.sin(osc.theta) + osc.gamma*osc.dtheta )


class BaseIntegrator:

    def __init__(self, _dt=0.01) :
        self.dt = _dt   # time step

    def integrate(self, simsystem, osc, obs):

        """ Perform a single integration step """
        
        self.timestep(simsystem, osc, obs)

        # Append observables to their lists
        obs.time.append(osc.t)
        obs.pos.append(osc.theta)
        obs.vel.append(osc.dtheta)
        # Function 'isinstance' is used to check if the instance of the system object is 'Harmonic' or 'Pendulum'
        if isinstance(simsystem, Harmonic):
            # Harmonic oscillator energy
            obs.energy.append(0.5 * osc.m * osc.L ** 2 * osc.dtheta ** 2 + 0.5 * osc.m * G * osc.L * osc.theta ** 2)
        else:
            # Pendulum energy
            obs.energy.append(0.5*osc.m*osc.L**2*osc.dtheta**2 + osc.m*G*osc.L*(1 - np.cos(osc.theta)))


    def timestep(self, simsystem, osc, obs):

        """ Virtual method: implemented by the child classes """
        
        pass


# HERE YOU ARE ASKED TO IMPLEMENT THE NUMERICAL TIME-MARCHING SCHEMES:

class EulerCromerIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel = simsystem.force(osc) / osc.m
        osc.t += self.dt
        osc.dtheta = osc.dtheta - (G/osc.L*np.sin(osc.theta) + osc.gamma*osc.dtheta)*self.dt
        osc.theta = osc.theta + osc.dtheta*self.dt
        


class VerletIntegrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        accel_n = simsystem.force(osc) / osc.m
        osc.t += self.dt
        osc.theta = osc.theta + osc.dtheta*self.dt + 1/2*accel_n*self.dt**2
        accel_n_p_1 = -(G/osc.L*np.sin(osc.theta) + osc.gamma*(osc.dtheta + self.dt/2*accel_n))/(1 + osc.gamma*self.dt/2)
        osc.dtheta = osc.dtheta + 1/2*(accel_n_p_1 + accel_n)*self.dt


class RK4Integrator(BaseIntegrator):
    def timestep(self, simsystem, osc, obs):
        osc.t += self.dt
        theta0 = osc.theta
        dtheta0 = osc.dtheta

        a1 = self.dt*simsystem.force(osc) / osc.m 
        b1 = self.dt*osc.dtheta

        osc.theta = theta0 + b1/2
        osc.dtheta = dtheta0 + a1/2
        a2 = self.dt*simsystem.force(osc)/osc.m 
        b2 = self.dt*osc.dtheta

        osc.theta = theta0 + b2/2
        osc.dtheta = dtheta0 + a2/2
        a3 = self.dt*simsystem.force(osc)/osc.m
        b3 = self.dt*osc.dtheta

        osc.theta = theta0 + b3
        osc.dtheta = dtheta0 + a3
        a4 = self.dt*simsystem.force(osc)/osc.m
        b4 = self.dt*osc.dtheta

        osc.theta = theta0 + 1/6*(b1 + 2*b2 + 2*b3 + b4)
        osc.dtheta = dtheta0 + 1/6*(a1 + 2*a2 + 2*a3 + a4)



# Animation function which integrates a few steps and return a line for the pendulum
def animate(framenr, simsystem, oscillator, obs, integrator, pendulum_line, stepsperframe):
    
    for it in range(stepsperframe):
        integrator.integrate(simsystem, oscillator, obs)

    x = np.array([0, np.sin(oscillator.theta)])
    y = np.array([0, -np.cos(oscillator.theta)])
    pendulum_line.set_data(x, y)
    return pendulum_line,


class Simulation:

    def reset(self, osc=Oscillator()) :
        self.oscillator = osc
        self.obs = Observables()

    def __init__(self, osc=Oscillator()) :
        self.reset(osc)

    # Run without displaying any animation (fast)
    def run(self,
            simsystem,
            integrator,
            tmax=30.,               # final time
            ):

        n = int(tmax / integrator.dt)

        for it in range(n):
            integrator.integrate(simsystem, self.oscillator, self.obs)

    # Run while displaying the animation of a pendulum swinging back and forth (slow-ish)
    # If too slow, try to increase stepsperframe
    def run_animate(self,
            simsystem,
            integrator,
            tmax=30.,               # final time
            stepsperframe=1         # how many integration steps between visualising frames
            ):

        numframes = int(tmax / (stepsperframe * integrator.dt))-2

        # WARNING! If you experience problems visualizing the animation try to comment/uncomment this line
        plt.clf()

        # If you experience problems visualizing the animation try to comment/uncomment this line
        # fig = plt.figure()

        ax = plt.subplot(xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
        plt.axhline(y=0)  # draw a default hline at y=1 that spans the xrange
        plt.axvline(x=0)  # draw a default vline at x=1 that spans the yrange
        pendulum_line, = ax.plot([], [], lw=5)
        plt.title(title)
        # Call the animator, blit=True means only re-draw parts that have changed
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[simsystem,self.oscillator,self.obs,integrator,pendulum_line,stepsperframe],
                                       frames=numframes, interval=25, blit=True, repeat=False)

        # If you experience problems visualizing the animation try to comment/uncomment this line
        # plt.show()

        # If you experience problems visualizing the animation try to comment/uncomment this line
        plt.waitforbuttonpress(30)

    # Plot coordinates and energies (to be called after running)
    def plot_observables(self, title="simulation", ref_E=None):

        #plt.clf()
        plt.title(title)
        plt.plot(self.obs.time, self.obs.pos, 'b-', label="Position")
        plt.plot(self.obs.time, self.obs.vel, 'r-', label="Velocity")
        plt.plot(self.obs.time, self.obs.energy, 'g-', label="Energy") # To see energy for large step sizes(self.obs.energy - np.ones(len(self.obs.energy))*4.3987893)*100000000
        if ref_E != None :
            plt.plot([self.obs.time[0],self.obs.time[-1]] , [ref_E, ref_E], 'k--', label="Ref.")
        plt.xlabel('time')
        plt.ylabel('observables')
        #plt.savefig(title + ".pdf")
        plt.legend()
        plt.show()

def plot_exact_harmonic(t_max, dt, init, w):
    n = int(t_max/dt)
    x = np.linspace(0, t_max, n+1)
    plt.title("Analytical harmonic solution")
    plt.plot(x, init*np.cos(w*x), "-b", label="Position")
    plt.plot(x, -init*w*sin(w*x), "-r", label="Velocity")
    plt.xlabel("time")
    plt.ylabel("observables")
    plt.legend()
    plt.show()

def plot_order(dt_list, t_max = 5):
    def get_RK4(t_max):
        S = Simulation()
        Integrator = RK4Integrator()
        Integrator.dt = 10**-6
        S.run(Pendulum(), Integrator, tmax=t_max)
        print(S.obs.time[-1])
        return S.obs.pos[-1]

    theta_diff_list = []
    theta_exact = get_RK4(t_max)
    S = Simulation()
    Integrator = VerletIntegrator()
    for dt in dt_list:
        S.reset(Oscillator())
        Integrator.dt = dt
        S.run(Pendulum(), Integrator, tmax=t_max)
        theta_temp = S.obs.pos[-1]
        print(S.obs.time[-1])
        theta_diff_list.append(abs(theta_temp - theta_exact))
    def get_p(E1, E2):
        return np.log(E1/E2)/np.log(2)
    x = dt_list[:-1]
    y = []
    for i in range(len(theta_diff_list)):
        y.append(get_p(theta_diff_list[i], theta_diff_list[i+1]))
        if i == len(theta_diff_list) - 2:
            break
    y = list(reversed(y))
    plt.scatter(x, y, label="p (time order of global error)")
    plt.axhline(y=2, color="r", linestyle="-", label="p=2")
    plt.ylim((0,5))
    plt.legend()
    plt.show()

def get_T(obs, tol=10**-3):
    dtheta = np.inf
    i = 0
    while True:
        try:
            if abs(obs.vel[i]) > tol:
                break
        except IndexError:
            return 0
        i += 1
    while abs(obs.vel[i]) > tol:
        i += 1
    return 2*obs.time[i]

def get_ODE(t, gamma, w):
    if gamma/2 > w:
        a = np.sqrt((gamma/2)**2 - w**2)
        return (1/2 + gamma/(4*a))*np.e**((-gamma/2 + a)*t) + (1/2 - gamma/(4*a))*np.e**((-gamma/2 - a)*t)
    elif gamma/2 == w:
        return (1 + gamma*t/2)*np.e**(-gamma/2*t)
    else:
        w_d = np.sqrt(w**2 - (gamma/2)**2)
        return np.e**(-gamma/2*t)*(np.cos(w_d*t) + gamma/(2*w_d)*np.sin(w_d*t))


def exercise_11() :
    S = Simulation() # set c = g/l = 9 and m = 1 in Oscillator 

    simsystem = Pendulum()
    integrator = VerletIntegrator()
    tmax = 25
    dt = 0.01
    c = S.oscillator.c
    init = S.oscillator.theta

    #S.run_animate(simsystem, integrator, tmax=tmax) 

    integrator.dt = dt
    #S.run(simsystem, integrator, tmax=tmax)
    #S.plot_observables()
    plt.show()
    #plot_exact_harmonic(tmax, dt, init, np.sqrt(c))

    plot_order([0.008, 0.004, 0.002, 0.001])

def exercise_12(): 
    t_max = 2
    dt = 0.0001

    S = Simulation()
    simsystem = Pendulum()
    integrator = VerletIntegrator()
    integrator.dt = dt
    c = S.oscillator.c

    inits = np.linspace(0.01, 0.5, 50)*np.pi
    peturbation_list = []
    def get_peturbation_T(init):
        return 2*np.pi/np.sqrt(c)*(1 + 1/16*init**2 + 11/3072*init**4 + 173/737280*init**6)
    for init in inits:
        peturbation_list.append(get_peturbation_T(init))

    T_list = []

    for init in inits:
        S.reset(Oscillator(theta0=init))
        S.run(simsystem, integrator, tmax=t_max)
        T_list.append(get_T(S.obs))
    plt.plot(inits, T_list, label="T(theta(0))")
    #plt.plot(inits, peturbation_list, label="Peturbation series T(theta(0))")
    plt.axhline(y=2*np.pi/np.sqrt(c), color = "r", label="T for harmonic oscillator")
    plt.xlabel("theta(0)")
    plt.ylabel("T for pendulum")
    plt.legend()
    plt.show()

def exercise_13():
    w = 3
    gamma = 3 # [0.5, 1, 2, 3]
    x0 = 1
    v0 = 0
    t_max = 15
    dt = 0.01

    S = Simulation()
    simsystem = Harmonic()
    integrator = VerletIntegrator()
    integrator.dt = dt

    x = np.linspace(0, t_max, int(t_max/dt)+10)
    #y = [get_ODE(t, gamma, w) for t in x]
    y = [np.e**(-gamma/2*t) for t in x]

    
    S.reset(Oscillator(c=w**2, gamma=gamma, theta0=x0, dtheta0=v0))
    S.run(simsystem, integrator, tmax=t_max)
    plt.plot(x, y, label="Envelope", color="C1")
    S.plot_observables()
    
        

if __name__ == "__main__" :
    #exercise_11()
    #exercise_12()
    exercise_13()
