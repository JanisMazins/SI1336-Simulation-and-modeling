import matplotlib.pyplot as plt
import random
import numpy as np
import os

def rw_points(steps=100, v0=(0,0), dx=1):
    x = v0[0]
    y = v0[1]
    rw_x = [x]
    rw_y = [y]
    rw_n = []
    for i in range(steps):
        n = random.randint(0,3)
        if n < 2:
            x0 = np.cos(np.pi*n)
            x += x0
        else:
            y0 = np.sin(np.pi/2 + np.pi*n)
            y += y0
        rw_x.append(x)
        rw_y.append(y)
        rw_n.append(n)
    return rw_x, rw_y, rw_n

def new_rw_points(steps=128, v0=(0,0), r0=1, a=1, c=1, m=100):
    x = v0[0]
    y = v0[1]
    rw_x = [x]
    rw_y = [y]
    rw_n = []
    for i in range(steps):
        r0 = (a*r0 + c) % m
        n = int(4*r0/m)
        if n < 2:
            x0 = np.cos(np.pi*n)
            x += x0
        else:
            y0 = np.sin(np.pi/2 + np.pi*n)
            y += y0
        rw_x.append(x)
        rw_y.append(y)
        rw_n.append(n)
    return rw_x, rw_y, rw_n
        
def plot_rw(rw_x, rw_y, rw_n, fig_name=None):
    plt.plot(rw_x[0], rw_y[0], "bx", label="starting point", markersize=10)
    plt.plot(rw_x[-1], rw_y[-1], "rx", label="end point", markersize=10)

    for i in range(len(rw_n)):
        shape_dict = {0:"right", 1:"left", 2:"right", 3:"left"}
        color_dict = {0:"r", 1:"aqua", 2:"b", 3:"black"}
        n = rw_n[i]
        plt.arrow(rw_x[i], rw_y[i], rw_x[i+1] - rw_x[i], rw_y[i+1] - rw_y[i], color=color_dict[n], shape=shape_dict[n], width=0.05)
        #plt.plot((rw_x[i],rw_x[i+1]), (rw_y[i], rw_y[i+1]), color=color_dict[n], linestyle=linestyle_dict[n])
    plt.figtext(0.7, 0.7, "red -> right", wrap=True, bbox={'facecolor':"red"})
    plt.figtext(0.7, 0.65, "teal -> up", wrap=True, bbox={'facecolor':"aqua"})
    plt.figtext(0.7, 0.6, "blue -> left", wrap=True, bbox={'facecolor':"blue"})
    plt.figtext(0.7, 0.55, "black -> down", wrap=True, color="white", bbox={'facecolor':"black"})
    plt.legend()
    plt.axis("equal")

    if fig_name is None:
        pass
    else:
        fig_path = r"C:\Users\Janis\coding\SI1336-Simulation-and-modeling\2_random_processes_and_complex_systems\figures" + "/" + fig_name + ".png"
        if os.path.exists(fig_path):
            new_fig_name = input(f"Proposed figure name {fig_name} already exists! Pick new name:")
            fig_path = r"C:\Users\Janis\coding\SI1336-Simulation-and-modeling\2_random_processes_and_complex_systems\figures" + "/" + new_fig_name + ".png"
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
    
    rw_x, rw_y, rw_n = rw_points(steps=steps, v0=v0)
    plot_rw(rw_x, rw_y, rw_n, fig_name=None)

def e_21b():
    steps = 100
    v0 = (0,0)
    r0 = 1
    a = 3
    c = 4
    m = 130

    rw_x, rw_y, rw_n = new_rw_points(steps=steps, v0=v0, r0=r0, a=a, c=c, m=m)
    plot_rw(rw_x, rw_y, rw_n, fig_name=None)#"new_rw_n100_r1_a3_c4_m128")


if __name__ == "__main__":
    #e_21a()
    e_21b()