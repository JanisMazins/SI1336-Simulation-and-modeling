import matplotlib.pyplot as plt
import random
import numpy as np
import os
import copy

def rw_points(get_n, args_n, steps=100, v0=(0,0), dx=1, mode=None):
    x = v0[0]
    y = v0[1]
    rw_x = [x]
    rw_y = [y]
    rw_n = []
    #Deals with special cases
    for i in range(steps):
        if get_n is get_recursive_n:
            r0, n = get_n(*args_n)
            args_n = (r0, *args_n[1:])
        elif get_n is random.choice:
            args_n = [[0, 1, 2, 3]]
            try:
                n_prev = rw_n[-1]
                match n_prev:
                    case 0:
                        args_n[0].remove(1)
                    case 1:
                        args_n[0].remove(0)
                    case 2:
                        args_n[0].remove(3)
                    case 3:
                        args_n[0].remove(2)
            except IndexError:
                pass

        n = get_n(*args_n)

        match n:
            case 0:
                x += 1
            case 1:
                x -= 1
            case 2:
                y += 1
            case 3:
                y -= 1

        if mode == "SA":
            for idx, xval in enumerate(rw_x):
                if x == xval and y == rw_y[idx]:
                    return None

        rw_x.append(x)
        rw_y.append(y)
        rw_n.append(n)
    return rw_x, rw_y, rw_n

def get_recursive_n(r0, a, c, m):
    r0 = (a*r0 + c) % m
    n = int(4*r0/m)
    return r0, n
        
def plot_rw(rw_x, rw_y, rw_n, fig_name=None):
    plt.plot(rw_x[0], rw_y[0], "bx", label="starting point", markersize=10)
    plt.plot(rw_x[-1], rw_y[-1], "rx", label="end point", markersize=10)

    for i in range(len(rw_n)):
        shape_dict = {0:"right", 1:"left", 2:"right", 3:"left"}
        color_dict = {0:"r", 1:"aqua", 2:"b", 3:"black"}
        n = rw_n[i]
        plt.arrow(rw_x[i], rw_y[i], rw_x[i+1] - rw_x[i], rw_y[i+1] - rw_y[i], color=color_dict[n], shape=shape_dict[n], width=0.05)
    plt.figtext(0.7, 0.7, "red -> right", wrap=True, bbox={'facecolor':"red"})
    plt.figtext(0.7, 0.65, "teal -> left", wrap=True, bbox={'facecolor':"aqua"})
    plt.figtext(0.7, 0.6, "blue -> up", wrap=True, bbox={'facecolor':"blue"})
    plt.figtext(0.7, 0.55, "black -> down", wrap=True, color="white", bbox={'facecolor':"black"})
    plt.legend()
    plt.axis("equal")

    if fig_name is None:
        pass
    else:
        current = os.getcwd()
        fig_path = current + r"\figures" + "/" + fig_name + ".png"
        if os.path.exists(fig_path):
            new_fig_name = input(f"Proposed figure name {fig_name} already exists! Pick new name:")
            if new_fig_name == "":
                raise("Run aborted.")
            fig_path = current + r"\figures" + "/" + new_fig_name + ".png"
        plt.savefig(fig_path, dpi=300)

    plt.show()

    if fig_name is not None:
        save = input("Save image? [y/n]")
        if save == "y":
            pass
        else:
            os.remove(fig_path)
        
def e_21a():
    steps = 100
    v0 = (0,0)
    
    rw_x, rw_y, rw_n = rw_points(random.randint, (0,3), steps=steps, v0=v0)
    plot_rw(rw_x, rw_y, rw_n, fig_name=None)

def e_21b():
    steps = 100
    v0 = (0,0)
    r0 = -100
    a = 3
    c = 4
    m = 128

    rw_x, rw_y, rw_n = rw_points(get_recursive_n, (r0, a, c, m), steps=steps, v0=v0)
    plot_rw(rw_x, rw_y, rw_n, fig_name="new_rw_n100_rm100_a3_c4_m128")

def e_21c():
    N_vals = 100
    N_max = 1000
    k = 1000

    N_list = np.round(np.linspace(1, N_max, N_vals)).astype(int)
    R2_rms_list = []
    std_list = []

    for N in N_list:
        print(f"N={N}/{N_max}")
        R2_temp_list = []
        R_temp_list = []
        for i in range(k):
            rw_x, rw_y, _ = rw_points(random.randint, (0,3), steps=N)
            R2 = rw_x[-1]**2 + rw_y[-1]**2
            R2_temp_list.append(R2)
            R_temp_list.append(np.sqrt(R2))
        mean = np.sqrt(sum(R2_temp_list)/k)
        var = sum(R2_temp_list)/k - (sum(R_temp_list)/k)**2
        std = np.sqrt(var/(k-1))
        R2_rms_list.append(mean)
        std_list.append(std)

    m1, b1 = np.polyfit(np.log(N_list.copy()), np.log(R2_rms_list.copy()), 1)
    f1 = lambda x: np.e**b1*x**m1

    plt.figure("Distance RMS")
    plt.plot(N_list, R2_rms_list, label="Root mean square distance")
    plt.plot(N_list, f1(N_list), label=f"y = e^[{round(b1, 2)}]*N^[{round(m1, 2)}]")
    plt.xlabel("N")
    plt.ylabel("Distance")
    plt.legend()

    m2, b2 = np.polyfit(np.log(N_list[1:].copy()), np.log(std_list[1:].copy()), 1)
    f2 = lambda x: np.e**b2*x**m2

    plt.figure("Distance STD")
    plt.plot(N_list, std_list, label="Standard deviation")
    plt.plot(N_list[1:], f2(N_list[1:]), label=f"y = e^[{round(b2, 2)}]*N^[{round(m2, 2)}]")
    plt.xlabel("N")
    plt.ylabel("Distance")
    plt.legend()

    plt.show()

def e_21d():
    N_vals = 30
    N_max = 30
    k = 1000000
    get_n_args = (0, 3)
    get_n = random.randint

    N_list = np.round(np.linspace(1, N_max, N_vals)).astype(int)
    rel_succ_list = []

    for N in N_list:
        print(f"N={N}/{N_max}")
        s = 0
        for i in range(k):
            try:
                _, _, _ = rw_points(get_n, get_n_args, steps=N, mode="SA")
                s += 1
            except TypeError:
                pass
        rel_succ_list.append(s/k)

    plt.figure("Success rate vs number of steps")
    plt.plot(N_list, rel_succ_list, label="Relative success vs nr. of steps N")
    plt.axhline(y=0, label="y=0", color="r", linewidth=0.4)
    plt.xlabel("N")
    plt.ylabel("Relative success")
    
    plt.legend()
    plt.show()

def e_21e():
    N_vals = 50
    N_max = 50
    k = 1000

    get_n = random.choice
    args_n = (0,3)
    mode = "SA"

    N_list = np.round(np.linspace(1, N_max, N_vals)).astype(int)
    R2_rms_list = []

    for N in N_list:
        print(f"N={N}/{N_max}")
        R2_temp_list = []
        for i in range(k):
            while True:
                try:
                    rw_x, rw_y, _ = rw_points(get_n, args_n, steps=N, mode=mode)
                    break
                except TypeError:
                    pass
            R2 = rw_x[-1]**2 + rw_y[-1]**2
            R2_temp_list.append(R2)
        mean = np.sqrt(sum(R2_temp_list)/k)
        R2_rms_list.append(mean)
    
    """
    A = np.log(N_list)
    b = np.log(R2_rms_list)

    m = solve(A.T@A, A.T@b)
    """
    m, b = np.polyfit(np.log(N_list.copy()), np.log(R2_rms_list.copy()), 1)
    f = lambda x: np.e**b*x**m

    plt.figure("Distance RMS")
    plt.plot(N_list, R2_rms_list, label="<R^2>^-(1/2)")
    plt.plot(N_list, f(N_list), label=f"y = e^[{round(b, 2)}]*N^[{round(m, 2)}]")
    plt.xlabel("N")
    plt.ylabel("Distance")
    plt.legend()

    plt.show()

def test_random():
    k = 100000
    n_list = [int(random.random()*4) for i in range(k)]
    plt.hist(n_list)
    plt.show()

def test_self_avoiding():
    steps = 150
    v0 = (0,0)
    
    i = 1

    while True:
        try:
            print(f"Iteration: {i}")
            i += 1
            rw_x, rw_y, rw_n = rw_points(random.choice, (0,3), steps=steps, v0=v0, mode="SA")
            break
        except TypeError:
            pass
    plot_rw(rw_x, rw_y, rw_n, fig_name=None)

if __name__ == "__main__":
    #e_21a()
    #e_21b()
    #e_21c()
    #e_21d()
    e_21e()
    #test_random()
    #test_self_avoiding()