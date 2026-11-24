import numpy as np
import matplotlib.pyplot as plt
import random

def m_step(x, delta, p):
    x0 = x
    d = delta*random.random()*random.choice([-1, 1])
    x1 = x + d

    r = random.random()
    if p(x1)/p(x0) <= r:
        x1 = x0

    return x1

def avg_f(points, f):
    s = [f(x) for x in points]
    N = len(points)
    return sum(s)/N

def metropolis(n_steps, delta, N_0=1000, initial_x=1, p=lambda x: x*np.e**(-x), f=lambda x: x):
    x = initial_x
    points = []

    for i in range(n_steps):
        x = m_step(x, delta, p)
        if i > N_0 - 1:
            points.append(x)

    avg = avg_f(points, f)
    
    return points, avg



def main():
    n_steps = 10**3
    deltas = np.logspace(-1, 1, 25)
    N_0 = 500
    k = 100

    rms_list = []
    std_list = []

    for delta in deltas:
        print(f"delta: {delta}/{deltas[-1]}")
        avg_list = []
        for i in range(k):
            _, avg = metropolis(n_steps, delta, N_0=N_0)
            avg_list.append(avg**2)
        rms = np.sqrt(sum(avg_list)/len(avg_list))
        std = np.sqrt(sum([(x - sum(avg_list)/len(avg_list))**2 for x in avg_list])/(len(avg_list) - 1))/np.sqrt(len(avg_list))
        rms_list.append(rms)
        std_list.append(std)
    
    plt.figure("RMS vs delta")
    plt.plot(deltas, rms_list, "r", label="RMS")
    plt.axhline(y=2, label="y=2")
    plt.xlabel("Delta")
    plt.ylabel("RMS")
    plt.legend()

    plt.figure("STD vs delta")
    plt.plot(deltas, std_list, "b", label="STD")
    plt.xlabel("Delta")
    plt.ylabel("STD")
    plt.legend()    

    plt.show()


if __name__ == "__main__":
    main()