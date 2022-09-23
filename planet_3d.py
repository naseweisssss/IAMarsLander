import numpy as np
import matplotlib.pyplot as plt

def Euler_3d(m, t_max, dt, r_init, v_init):
    # initialise a time array
    t_array = np.arange(0, t_max, dt)

    # initialise empty lists to record trajectories
    r_list = []
    v_list = []

    r = r_init
    v = v_init
    for t in t_array:
        r_list.append(r)
        v_list.append(v)

        a = Force(m,r)/m
        r = r + dt*v
        v = v + dt*a
    
    return [t_array, r_list, v_list]

def Verlet_3d(m, t_max, dt, r_init, v_init):
    # initialise a time array
    t_array = np.arange(0, t_max, dt)

    # initialise empty lists to record trajectories
    r_list = []
    v_list = []

    r = r_init
    v = v_init
    for t in t_array:
        r_list.append(r)
        v_list.append(v)

        r_past = r

        try:
            r = 2*r - r_list[-2] + (dt**2)*Force(m, r_list[-1])/m
        except:
            r = r + dt*v
        v = (r - r_past)/dt
    
    return [t_array, r_list, v_list]

def Force(m, r):
    G = 6.6742e-11
    M = 6.42e23
    r_mag = np.linalg.norm(r)
    r_hat = r/r_mag
    F = -(G*M*m)/(r_mag**2)*r_hat
    return F

def simulate(m, t_max, dt, r_init, v_init, plot_type):
    t_array, r_list, v_list = Euler_3d(m, t_max, dt, r_init, v_init)
    t_array_vel, r_list_vel, v_list_vel = Verlet_3d(m, t_max, dt, r_init, v_init)
    if plot_type == "altitude":
        fig = plt.figure(1)
        ax1 = plt.subplot(211)
        ax1.set_title(" Euler Integration Results")
        ax1.set_xlabel('time (s)')
        ax1.grid()
        ax1.plot(t_array, [pos[2] for pos in r_list], label='x (m)')
        # ax1.plot(t_array, v_list, label='v (m/s)')
        ax1.legend()
        ax2 = plt.subplot(212)
        ax2.set_title("Verlet Integration Results")
        ax2.set_xlabel('time (s)')
        ax2.grid()
        ax2.plot(t_array_vel,  [pos[2] for pos in r_list_vel], label='x (m)')
        # ax2.plot(t_array_vel, v_list_vel, label='v (m/s)')
        ax2.legend()
        fig.tight_layout()
        plt.show()
    elif plot_type == "trajectory":
        print(r_list)
        fig = plt.figure(1)
        ax1 = plt.subplot(211)
        ax1.set_title(" Euler Integration Results")
        ax1.set_xlabel('y_coordinate (s)')
        ax1.grid()
        ax1.plot([pos[1] for pos in r_list], [pos[2] for pos in r_list], label='x (m)')
        ax1.plot(0,0, "ro", label = "origin")
        # ax1.plot(t_array, v_list, label='v (m/s)')
        ax1.legend()
        ax2 = plt.subplot(212)
        ax2.set_title("Verlet Integration Results")
        ax2.set_xlabel('y_coordinate (s)')
        ax2.grid()
        ax2.plot([pos[1] for pos in r_list_vel],  [pos[2] for pos in r_list_vel], label='x (m)')
        ax2.plot(0,0, "ro", label = "origin")
        # ax2.plot(t_array_vel, v_list_vel, label='v (m/s)')
        ax2.legend()
        fig.tight_layout()
        plt.show()

if __name__ == "__main__":
    # Simulate straight down descent
    simulate(m=10, t_max = 10, dt =  0.1, 
    r_init = np.array([0,0,3.3e6], dtype= np.float64),
    v_init = np.array([0,0,0], dtype = np.float64),
    plot_type = "altitude")
    #Simulate Elliptical orbit
    simulate(m=10, t_max = 2050, dt =  1, 
    r_init = np.array([0,0,3.3e6], dtype= np.float64),
    v_init = np.array([0,1000,0], dtype = np.float64),
    plot_type = "trajectory")
    #Simuate Circular Orbit
    G = 6.6742e-11
    M = 6.42e23
    R = 3.3e6
    orbit_vel = np.sqrt(G*M/R)
    simulate(m=10, t_max = 7036, dt = 2, 
    r_init = np.array([0,0,R], dtype= np.float64),
    v_init = np.array([0,orbit_vel,0], dtype = np.float64),
    plot_type = "trajectory")
    #Simulate Hyperbolic escape
    escape_vel = orbit_vel*np.sqrt(2)
    simulate(m=10, t_max = 10000, dt =  2, 
    r_init = np.array([0,0,3.3e6], dtype= np.float64),
    v_init = np.array([0,escape_vel + 100,0], dtype = np.float64),
    plot_type = "trajectory")

