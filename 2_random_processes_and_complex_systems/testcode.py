def e_22c():
    vmax = 2
    p = 0.5
    max_roadLength = 75
    min_roadLength = 2

    S = Simulation()
    
    N_stride = 1
    k = 100
    numsteps = 100

    fr_max_list = []
    fr_max_idx_list = []

    for roadLength in range(min_roadLength, max_roadLength + 1):
        print(f"road length = {roadLength}/{max_roadLength}")
        fr_list = []

        for N_cars in range(0, roadLength, N_stride):
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
        fr_max_list.append(fr_max)
        fr_max_idx_list.append(fr_max_idx/roadLength)

    roadLengths = [i for i in range(min_roadLength, max_roadLength + 1)]

    plt.plot(roadLengths, fr_max_list, label="Maximum flow rate", color="b")
    plt.plot(roadLengths, fr_max_idx_list, label="Density at maximum flow rate", color="r")
    plt.legend()
    plt.xlabel(f"Roadlength")
    plt.ylabel("Observables")
    plt.show()