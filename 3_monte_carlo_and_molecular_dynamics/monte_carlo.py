import numpy as np
import matplotlib.pyplot as plt
import random

def m_step(x, delta, p, tracking_success_rate=False):
    x0 = x
    d = delta*random.random()*random.choice([-1, 1])
    x1 = x + d
    success = True

    r = random.random()
    if p(x1)/p(x0) <= r:
        x1 = x0
        success = False

    if tracking_success_rate:
        return x1, success

    return x1

def avg_f(points, f):
    s = [f(x) for x in points]
    N = len(points)
    return sum(s)/N

def metropolis(n_steps, delta, N_0=100, initial_x=1, p=lambda x: x*np.e**(-x), f=lambda x: x, tracking_success_rate=False):
    x = initial_x
    points = []
    successes = 0

    for i in range(n_steps + N_0):
        if tracking_success_rate: #Tracking success rate
            x, success = m_step(x, delta, p, tracking_success_rate=tracking_success_rate)
        else: #Not tracking success rate
            x = m_step(x, delta, p)
            success = False

        if i > N_0 - 1: #Only take points after N_0 
            points.append(x)
        
        if success:
            successes += 1

    avg = avg_f(points, f)
    
    if tracking_success_rate:
        sr = successes/(n_steps + N_0)
        return points, avg, sr
    else:
        return points, avg



def rms_std():
    n_steps = 10**4
    N_0 = 1000
    p = lambda x: x*np.e**(-x)

    deltas_log = np.logspace(-1, 2, 10)
    deltas_lin = np.linspace(0.1, 1e2, 10)
    deltas = np.unique(np.sort(np.concatenate((deltas_log, deltas_lin))))
    
    k = 100

    rms_list = []
    std_list = []

    for delta in deltas:
        print(f"delta: {delta}/{deltas[-1]}")
        avg_list = []
        for i in range(k):
            _, avg = metropolis(n_steps, delta, N_0=N_0)
            avg_list.append((avg-2)**2)
        rms = np.sqrt(sum(avg_list)/len(avg_list))
        std = np.sqrt(sum([(x - sum(avg_list)/len(avg_list))**2 for x in avg_list])/(len(avg_list) - 1))/np.sqrt(len(avg_list))
        rms_list.append(rms)
        std_list.append(std)
    
    plt.figure("RMS vs delta")
    plt.plot(deltas, rms_list, "r", label="RMS")
    plt.axhline(y=0, label="y=0")
    plt.xlabel("Delta")
    plt.ylabel("RMS")
    plt.legend()

    plt.figure("STD vs delta")
    plt.loglog(deltas, std_list, "b", label="STD")
    plt.xlabel("Delta")
    plt.ylabel("STD")
    plt.legend()    

    x = np.linspace(0, 15, 10**5)

    plt.figure("x*e^-x")
    plt.plot(x, p(x))
    plt.show()

def success_rate():
    deltas_log = np.logspace(-1, 2, 10)
    deltas_lin = np.linspace(0.1, 10**2, 10)
    deltas = np.unique(np.sort(np.concatenate((deltas_log, deltas_lin))))

    n_steps = 10**5
    N_0 = 1000
    p = lambda x: x*np.e**(-x)

    sr_list = []

    for delta in deltas:
        print(f"delta: {delta}/{deltas[-1]}")
        _, _, sr = metropolis(n_steps, delta, N_0=N_0, tracking_success_rate=True)
        sr_list.append(sr)

    m, b = np.polyfit(np.log([d for d in deltas if d >= 4.5]), np.log([sr for idx, sr in enumerate(sr_list) if deltas[idx] > 4.5]), 1)
    g = lambda x: np.e**b*x**m

    plt.figure("Success rate of moves")
    plt.loglog(deltas, sr_list, label="Success rate")
    plt.plot(deltas, g(deltas), label=f"y = e^[{round(b, 2)}]*delta^[{round(m, 2)}]")
    plt.xlabel("Delta")
    plt.ylabel("Success rate")
    plt.legend()
    plt.show()



if __name__ == "__main__":
    rms_std()
    #success_rate()