import matplotlib.pyplot as plt
import random
import numpy as np
import os

def rw_points(get_n, args_n, steps=100, v0=(0,0), dx=1):
    x = v0[0]
    y = v0[1]
    rw_x = [x]
    rw_y = [y]
    rw_n = []
    for i in range(steps):
        if get_n is get_recursive_n:
            r0, n = get_n(*args_n)
            args_n = (r0, *args_n[1:])
        else:
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
    plt.figtext(0.7, 0.65, "teal -> up", wrap=True, bbox={'facecolor':"aqua"})
    plt.figtext(0.7, 0.6, "blue -> left", wrap=True, bbox={'facecolor':"blue"})
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
    N_vals = 500
    N_max = 2000
    k = 2000

    N_list = np.round(np.linspace(1, N_max, N_vals)).astype(int)
    R_mean_list = []

    for N in N_list:
        print(f"N={N}/{N_max}")
        R2_temp_list = []
        for i in range(k):
            rw_x, rw_y, _ = rw_points(random.randint, (0,3), steps=N)
            R2 = rw_x[-1]**2 + rw_y[-1]**2
            R2_temp_list.append(R2)
        mean = np.sqrt(sum(R2_temp_list)/k)
        R_mean_list.append(mean)

    plt.plot(N_list, R_mean_list, label="<R^2>^-(1/2)[N]")
    plt.xlabel("N")
    plt.ylabel("<R^2>^-(1/2)")
    plt.legend()
    plt.show()

    



if __name__ == "__main__":
    #e_21a()
    #e_21b()
    e_21c()