import pylab as plt
import numpy as np


class simulation():

    def __init__(self):
        self.R=8.31446261815324 #gas constant
        self.mu = 1 #TODO : find value of mu (nils)
        self.P_atm = 100000
        self.mu_G = 1
        self.eps = 5/2 #degrees of freedom gas
        self.T_ob = 300 #temperature outside
        self.V_p = 0.01
        self.L = 0.2 #legth of pipe
        self.molemass = 0.029 #mass of 1 mole of air
        self.airdensity = 1.29 #kg/m3
        self.a=2/3 #gas konstant for air

        self.t_list=[]
        self.T_list = []
        self.P_list = []
        self.V_list = []
        self.N_list = []

        self.P_prime=[]
        self.T_prime = []
        self.V_prime = []
        self.N_prime = []

        self.W=0

        # TODO: change parameters and see how results change
        # constants
        self.time_cylinder = 0.1  # time during which culinder moves (play around to see different results, less time - faster piston)
        self.cylinderheight = 0.15  # more distance - faster piston
        self.cylinderArea = 0.07 ** 2  # bigger area, bigger volume, more W
        self.V0 = self.cylinderArea * 0.03  # initial volume (play around to see different results)
        self.holeArea = 0.03 ** 2  # area of inlet (play around to see different results)
        self.t_tot = self.time_cylinder + 0.1  # total time simulation runs
        self.dt=0.00001


    endTime = 1

    #returns derivative of volume
    #piston has position (1-cos(pi*t/t_cylinder)*cylinderheight (tigonmetric)
    #one can easily change to a linear piston
    def volume_prime(self, time, cylinderArea, t_cylinder, cylinderheight):
        if t_cylinder<=time:
            return 0

        v_prime=np.sin(np.pi*(time/t_cylinder))*(cylinderheight/2)*cylinderArea/t_cylinder*np.pi
        return v_prime


    def pressure_prime(self, N_prime, P, T, V, V_prime, T_prime, N):
        '''Returning the analytical time derivative of the pressure, derived from the ideal gas law.'''
        #return ((N_prime * R * T) / V) - P  * (V_prime / V - T_prime / T)
        return -self.R*N*T*V_prime/(V**2) + self.R*(N*T_prime+N_prime*T)/V    #alternate way of comouting

    def N_prime_hole(self, P1, P2, T1, T2, m, holeArea):
        '''Georgios' Solution of particle flow from hole approximation'''
        return holeArea /( np.sqrt(12 * self.R * m) )* (P1 / np.sqrt(T1) - P2 / np.sqrt(T2))

    def N_prime_pipe(self, L, P, holeArea):
        R=np.sqrt(holeArea/np.pi)
        '''Nils / Masato's solution from cylindrical pipe approximation'''
        return -(1/self.V_p) * (np.pi / 8 * self.mu) * R**4 * ((P - self.P_atm) / L)

    def N_prime_simplepipe(self, P, holearea):
        "max solution for flow trough pipe without friction"
        if P>self.P_atm:
            P=self.P_atm
        u=np.sqrt(2*(self.P_atm-P)/self.airdensity)
        return u*holearea*self.airdensity/self.molemass

    def temp_prime(self, V_prim, P, T, N_prim, N):
        '''Returning the temperature dervied from the First Law of Thermodynamics (?)'''
        return (V_prim*(P-self.P_atm))/(self.a*self.R*N) + N_prim*(T)/N

    def temp_prime2(self, V, V_prim, P, P_prim, N, N_prim, T):
        '''Returning the temperature dervied from the First Law of Thermodynamics (?)'''
        return N_prim*(self.T_ob-T)/N - (V_prim*P)/(self.a*self.R)

    #returns numver of particles
    def N_func(self, p, v, T):
        return (p*v)/(T * self.R)

    def Euler_method(self):
        '''Performing the Euler Method on P, T, V and the work, W'''

        nsteps = int(self.t_tot / self.dt)


        #boundary coditions
        self.t_list=[0]
        self.P_list= [self.P_atm]
        self.T_list = [300]
        self.V_list = [self.V0]
        self.N_list = [self.N_func(self.P_list[0],self.V_list[0], self.T_list[0])]

        #boundary dynamic condition
        self.N_prime = [0]
        self.V_prime = [self.volume_prime(self.t_list[0], self.cylinderArea, self.time_cylinder, self.cylinderheight)]
        self.P_prime = [0]
        self.T_prime = [0]


        for i in range(1, nsteps):

            self.t_list.append(i*self.dt)

            # TODO: change pipes and see which work and if so, how they affect result
            #N_prime.append(N_prime_simplepipe(P_list[i - 1], holeArea))
            #N_prime.append(N_prime_hole(P_list[i-1], P_atm, T_list[i-1], T_ob, molemass, holeArea))
            self.N_prime.append(self.N_prime_pipe(self.L, self.P_list[i-1], self.holeArea))

            self.V_prime.append(self.volume_prime(self.t_list[i], self.cylinderArea, self.time_cylinder, self.cylinderheight))

            #T_prime.append(temp_prime(V_prime[-1], P_list[-1], T_list[-1], N_prime[-1], N_list[-1]))
            self.T_prime.append(self.temp_prime2(self.V_list[i-1], self.V_prime[i-1], self.P_list[i-1], self.P_prime[i-1], self.N_list[i-1], self.N_prime[i-1], self.T_list[i-1]))
            #T_prime.append(0)

            self.P_prime.append(self.pressure_prime(self.N_prime[i-1], self.P_list[i-1], self.T_list[i-1], self.V_list[i-1], self.V_prime[i-1], self.T_prime[i-1], self.N_list[i-1]))

            #euler step
            Ptemp=self.P_list[i-1]+self.dt*self.P_prime[i]
            self.P_list.append(Ptemp)

            Ttemp=self.T_list[i-1]+ self.dt*self.T_prime[i]
            self.T_list.append(Ttemp)

            Ntemp = self.N_list[i - 1] + self.dt * self.N_prime[i]
            self.N_list.append(Ntemp)

            Vtemp = self.V_list[i - 1] + self.dt * self.V_prime[i]
            self.V_list.append(Vtemp)

        for i in range(len(self.t_list)):
            self.W += self.V_prime[i] * (self.P_atm - self.P_list[i]) * self.dt


    def plotPresure(self):
        #ploting result
        plt.plot(self.t_list, self.P_list)
        plt.ylim([min(self.P_list)-1000, 101000])
        plt.scatter(self.time_cylinder, self.P_list[int(self.time_cylinder/self.dt)], label="time at which piston has reached bottom position")
        plt.legend()
        plt.xlabel("time")
        plt.ylabel("pressure (pascal)")
        plt.show()


    def plotTemp(self):
        plt.plot(self.t_list, self.T_list)
        plt.ylim([min(self.T_list) - 10, max(self.T_list)+10])
        plt.xlabel("time")
        plt.ylabel("temperature (kelvin)")
        plt.scatter(self.time_cylinder, self.T_list[int(self.time_cylinder / self.dt)],
                    label="time at which piston has reached bottom position")
        plt.legend()
        plt.show()

    def plotVolume(self):
        plt.plot(self.t_list, self.V_list)
        plt.xlabel("time")
        plt.ylabel("Volume (m3)")
        plt.scatter(self.time_cylinder, self.V_list[int(self.time_cylinder / self.dt)],
                    label="time at which piston has reached bottom position")
        plt.legend()
        plt.show()

    def plotN(self):
        plt.plot(self.t_list, self.N_list)
        plt.scatter(self.time_cylinder, self.N_list[int(self.time_cylinder / self.dt)],
                    label="time at which piston has reached bottom position")
        plt.legend()
        plt.xlabel("time")
        plt.ylabel("number of particles in moles")
        plt.show()


sim=simulation()
sim.Euler_method()
sim.plotPresure()
print(sim.W)


