#!/bin/python3

# Template for traffic simulation
# BH, MP 2021-11-15, latest version 2024-11-08.

"""
    This template is used as backbone for the traffic simulations.
    Its structure resembles the one of the pendulum project, that is you have:
    (a) a class containing the state of the system and it's parameters
    (b) a class storing the observables that you want then to plot
    (c) a class that propagates the state in time (which in this case is discrete), and
    (d) a class that encapsulates the aforementioned ones and performs the actual simulation
    You are asked to implement the propagation rule(s) corresponding to the traffic model(s) of the project.
"""

import math
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy.random as rng
import numpy as np
import random
import copy

import matplotlib


class Cars:

    """ Class for the state of a number of cars """

    def __init__(self, numCars=5, roadLength=50, v0=1, init_cond="basic"):
        self.numCars    = numCars
        self.roadLength = roadLength
        self.t  = 0
        self.x  = []
        self.v  = []
        self.c  = []
        positions = [i for i in range(self.roadLength)]
        for i in range(self.numCars):
            if init_cond == "basic":
                self.x.append(i)
            elif init_cond == "random":
                random_idx = random.randint(0, len(positions) - 1)
                random_x = positions.pop(random_idx)
                self.x.append(random_x)
            elif init_cond == "even":
                self.x.append((i*(self.roadLength/self.numCars)) % self.roadLength)
            self.v.append(v0)       # the speed of the cars
            self.c.append(i)        # the color of the cars (for drawing)
        self.x.sort()

    def distance(self, i):
        if i == self.numCars - 1:
            next_idx = 0
        else:
            next_idx = i + 1
        
        if self.x[next_idx] - self.x[i] > 0:
            return self.x[next_idx] - self.x[i]
        else:
            return self.roadLength - self.x[i] + self.x[next_idx]

    def fuel(self, i):
        v = self.v[i]
        return 1 + v + 0.2*v**3

class Observables:

    """ Class for storing observables """

    def __init__(self):
        self.time = []          # list to store time
        self.flowrate = []      # list to store the flow rate
        self.f = []
        

class BasePropagator:

    def __init__(self):
        return
        
    def propagate(self, cars, obs):

        """ Perform a single integration step """
        
        fr = self.timestep(cars, obs)
        try:
            avg_fuel = sum([cars.fuel(i) for i in range(cars.numCars)])/cars.numCars
        except ZeroDivisionError:
            avg_fuel = 0
        # Append observables to their lists
        obs.time.append(cars.t)
        obs.flowrate.append(fr)
        obs.f.append(avg_fuel)
         
              
    def timestep(self, cars, obs):

        """ Virtual method: implemented by the child classes """
        
        pass
      
        
class ConstantPropagator(BasePropagator) :
    
    """ 
        Cars do not interact: each position is just 
        updated using the corresponding velocity 
    """
    
    def timestep(self, cars, obs):
        for i in range(cars.numCars):
            cars.x[i] += cars.v[i]
        cars.t += 1
        return 0

class MyPropagator(BasePropagator) :

    def __init__(self, vmax, p):
        BasePropagator.__init__(self)
        self.vmax = vmax
        self.p = p
        p_bools = [True for i in range(round(1000*self.p))]
        p_complementary_bools = [False for i in range(1000 - round(1000*self.p))]
        self.bool_list = p_bools + p_complementary_bools

    def timestep(self, cars, obs):
        for i in range(cars.numCars):
            if cars.v[i] < self.vmax:
                cars.v[i] += 1 

            d = cars.distance(i)
            if cars.v[i] >= d:
                cars.v[i] = d - 1

            randomly_slower = random.choice(self.bool_list)
            if randomly_slower and cars.v[i] > 0:
                cars.v[i] -= 1

            if cars.x[i] + cars.v[i] < cars.roadLength:
                cars.x[i] += cars.v[i]
            else:
                cars.x[i] = cars.x[i] + cars.v[i] - cars.roadLength

        cars.t += 1
        fr = sum(cars.v)/cars.roadLength

        return fr

############################################################################################

def draw_cars(cars, cars_drawing):

    """ Used later on to generate the animation """
    theta = []
    r     = []

    for position in cars.x:
        # Convert to radians for plotting  only (do not use radians for the simulation!)
        theta.append(position * 2 * math.pi / cars.roadLength)
        r.append(1)

    return cars_drawing.scatter(theta, r, c=cars.c, cmap='hsv')


def animate(framenr, cars, obs, propagator, road_drawing, stepsperframe):

    """ Animation function which integrates a few steps and return a drawing """

    for it in range(stepsperframe):
        propagator.propagate(cars, obs)

    return draw_cars(cars, road_drawing),


class Simulation:

    def reset(self, cars=Cars()) :
        self.cars = cars
        self.obs = Observables()

    def __init__(self, cars=Cars()) :
        self.reset(cars)

    def plot_observables(self, title="simulation"):
        plt.clf()
        plt.title(title)
        plt.plot(self.obs.time, self.obs.flowrate)
        plt.xlabel('time')
        plt.ylabel('flow rate')
        #plt.savefig(title + ".pdf")
        plt.figure("Flow rate histogram")
        plt.hist(self.obs.flowrate, bins=int(np.sqrt(len(self.obs.flowrate))))
        plt.xlabel("flow rate")
        plt.show()

    # Run without displaying any animation (fast)
    def run(self,
            propagator,
            numsteps=200,           # final time
            title="simulation",     # Name of output file and title shown at the top
            ):

        for it in range(numsteps):
            propagator.propagate(self.cars, self.obs)

    # Run while displaying the animation of bunch of cars going in circe (slow-ish)
    def run_animate(self,
            propagator,
            numsteps=200,           # Final time
            stepsperframe=1,        # How many integration steps between visualising frames
            title="simulation",     # Name of output file and title shown at the top
            ):

        numframes = int(numsteps / stepsperframe)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='polar')
        ax.axis('off')
        # Call the animator, blit=False means re-draw everything
        anim = animation.FuncAnimation(plt.gcf(), animate,  # init_func=init,
                                       fargs=[self.cars,self.obs,propagator,ax,stepsperframe],
                                       frames=numframes, interval=50, blit=True, repeat=False)
        plt.show()

        # If you experience problems visualizing the animation and/or
        # the following figures comment out the next line 
        #plt.waitforbuttonpress(30)
    

# It's good practice to encapsulate the script execution in 
# a main() function (e.g. for profiling reasons)
def main() :
    vmax = 2
    p = 0.5
    numcars = 25
    roadlength = 50
    numsteps = 5000

    # Be sure you are passing the correct initial conditions!
    cars = Cars(numCars = numcars, roadLength=roadlength, init_cond="even")

    # Create the simulation object for your cars instance:
    S = Simulation(cars)

    #S.run_animate(propagator=ConstantPropagator())
    #S.run_animate(propagator=MyPropagator(vmax=vmax, p=p))

    S.run(propagator=MyPropagator(vmax=vmax, p=p), numsteps=numsteps)
    S.plot_observables()

def e_22a():
    vmax = 2
    p = 0.5
    roadLength = 200

    S = Simulation()
    
    max_N_cars = 200
    N_stride = 1
    k = 10
    numsteps = 100

    fr_list = []

    for N_cars in range(0, max_N_cars, N_stride):
        print(f"N_cars = {N_cars}/{max_N_cars}")
        fr_temp_list = []
        for i in range(k):
            S.reset(Cars(numCars=N_cars, roadLength=roadLength, init_cond="random"))
            S.run(propagator=MyPropagator(vmax=vmax, p=p), numsteps=numsteps)
            avg_fr_init_cond = sum(S.obs.flowrate)/S.cars.t
            fr_temp_list.append(avg_fr_init_cond)
        avg_fr_k_attempts = sum(fr_temp_list)/k
        fr_list.append(avg_fr_k_attempts)

    fr_max = 0
    fr_max_idx = None
    for idx, fr in enumerate(fr_list):
        if fr > fr_max:
            fr_max = fr
            fr_max_idx = idx
    N_list = [i for i in range(0, max_N_cars, N_stride)]
    densities = [N/roadLength for N in N_list]

    plt.plot(densities, fr_list, label="Flowrate", color="b")
    plt.axvline(x=densities[fr_max_idx], label=f"x={densities[fr_max_idx]} (fr_max = {fr_max})", color="r")
    plt.legend()
    plt.xlabel(f"Car density (number of cars / road length, road length = {roadLength}")
    plt.ylabel("Flow rate")
    plt.show()

def e_22b():
    vmax = 2
    p = 0.5
    roadLength = 50
    N_cars = 25

    S = Simulation()

    k = 10000
    numsteps = 100
    tol = 0.001
    fr_list = []
    fr2_list = []
    std_list = []
    last_idx = k - 1

    for i in range(k):
        S.reset(Cars(numCars=N_cars, roadLength=roadLength, init_cond="random"))
        S.run(propagator=MyPropagator(vmax=vmax, p=p), numsteps=numsteps)
        avg_fr = sum(S.obs.flowrate)/S.cars.t
        avg_fr2 = sum([fr**2 for fr in S.obs.flowrate])/S.cars.t
        fr_list.append(avg_fr)
        fr2_list.append(avg_fr2)

        try:
            var = sum(fr2_list)/len(fr2_list) - (sum(fr_list)/len(fr_list))**2
            std = np.sqrt(var/(len(fr_list) - 1)) #Standard Error Estimate
            std_list.append(std)
            print(f"iteration: {i+1}, std={std}")
            if std < tol and i > 0:
                last_idx = i
                break
        except ZeroDivisionError:
            std_list.append(0)
            print(f"iteration: {i+1}, std={0}")


    k_list = [j for j in range(last_idx + 1)]

    m, b = np.polyfit(np.log(k_list[1:].copy()), np.log(std_list[1:].copy()), 1)
    f = lambda x: np.e**b*x**m

    print(f"average flow rate: {sum(fr_list)/len(fr_list)}")
    
    plt.plot(k_list[1:], std_list[1:], label=f"Standard deviation per iteration ({last_idx + 1} iterations, mean = {round(sum(fr_list)/len(fr_list), 3)})", color="b")
    plt.plot(k_list[1:], f(k_list[1:]), label=f"y = [{round(np.e**b, 3)}*N^[{round(m, 3)}]")
    plt.legend()
    plt.xlabel(f"Number of attempts")
    plt.ylabel("Standard deviation")

    plt.figure()
    plt.hist(fr_list, bins=int(np.sqrt(len(fr_list))))

    plt.show()
    
def e_22c():
    vmax = 2
    p = 0.5

    min_roadLength = 2
    max_roadLength = 100
    N_stride = 1
    densities = np.linspace(0.1, 0.9, 5)

    S = Simulation()

    k = 10
    numsteps = 1000
    tol = 0.0001
    
    nested_fr_list = []
    diff_list = []

    for rL_idx, roadLength in enumerate(range(min_roadLength, max_roadLength + 1, N_stride)):
        try:
            print(f"road length = {roadLength}/{max_roadLength}, last diff = {diff_list[rL_idx - 1]}")
        except IndexError:
            print(f"road length = {roadLength}/{max_roadLength}")

        fr_list = []
        for d in densities:
            fr_temp_list = []
            for i in range(k):
                S.reset(Cars(numCars=int(roadLength*d), roadLength=roadLength, init_cond="random"))
                S.run(propagator=MyPropagator(vmax=vmax, p=p), numsteps=numsteps)
                avg_fr = sum(S.obs.flowrate)/S.cars.t
                fr_temp_list.append(avg_fr)
            fr_list.append(sum(fr_temp_list)/k)
        nested_fr_list.append(fr_list)
        try:
            diffs = [abs(nested_fr_list[rL_idx][j] - nested_fr_list[rL_idx - 1][j]) for j in range(len(densities))]
        except IndexError:
            diffs = [np.inf for i in range(k)]
        
        mean_diff = sum(diffs)/len(diffs)

        diff_list.append(mean_diff)

        if mean_diff < tol and rL_idx != 0 and True == False:
            max_roadLength = roadLength
            break

    plt.loglog(list(range(min_roadLength + N_stride, max_roadLength + 1, N_stride)), diff_list[1:], label="Average difference between fundemental diagram iterations")
    plt.xlabel("Road length")
    plt.ylabel("Mean difference between fundemental diagram iterations")
    plt.legend()

    plt.figure()
    roadLengths = list(range(min_roadLength, max_roadLength + 1, N_stride))
    for idx, d in enumerate(densities):
        plt.plot(roadLengths, [nested_fr_list[i][idx] for i in range(len(roadLengths))], label=f"d={round(d, 3)}")
        plt.legend()
    plt.xlabel("Road length")
    plt.ylabel("Flow rate")
    plt.legend()
    plt.show()

def e_22d():
    vmax_array = np.linspace(1, 50, 50)
    p = 0.25
    d = 0.4
    roadLength = 250
    k = 100
    
    numsteps = 1000

    fr_list = []
    fuel_list = []

    S = Simulation()

    for vmax in vmax_array:
        fr_temp_list = []
        fuel_temp_list = []
        for i in range(k):
            S.reset(Cars(numCars=int(d*roadLength), roadLength=roadLength, init_cond="random"))
            S.run(propagator=MyPropagator(vmax=vmax, p=p), numsteps=numsteps)

            avg_fr  = sum(S.obs.flowrate)/S.cars.t
            avg_fuel = sum(S.obs.f)/S.cars.t

            fr_temp_list.append(avg_fr)
            fuel_temp_list.append(avg_fuel)

        avg_avg_fr = sum(fr_temp_list)/k
        avg_avg_fuel = sum(fuel_temp_list)/k

        fr_list.append(avg_avg_fr)
        fuel_list.append(avg_avg_fuel)
    
    fuel_div_fr_list = [fuel_list[i]/fr_list[i] for i in range(len(vmax_array))]

    plt.figure()
    plt.plot(vmax_array, fr_list, color="r", label="Flow rate")
    plt.xlabel("Maximum speed")
    plt.ylabel("Flow rate")
    plt.legend()

    plt.figure()
    plt.plot(vmax_array, fuel_div_fr_list, color="b", label="Fuel consumption divided by flow rate")
    plt.xlabel("Maximum speed")
    plt.ylabel("Average fuel consumption / flow rate")
    plt.legend()

    plt.show()

# Calling 'main()' if the script is executed.
# If the script is instead just imported, main is not called (this can be useful if you want to
# write another script importing and utilizing the functions and classes defined in this one)
if __name__ == "__main__" :
    #main()
    #e_22a()
    #e_22b()
    e_22c()
    #e_22d()