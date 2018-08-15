"""
Run by:
>module load simpy
>python3 SimPysample.py

Basic program that gives you a brief overview of how
SimPy works.

"""

import simpy


#Simpy Sophia
#-----------------------------------------------------------
print("-->Sophias Beispiele: \n")
class Car(object):
    def __init__(self, env):
        self.env = env
        #starts the process everythime, when instance is created:
        self.action = env.process(self.runner()) 
        
    def runner(self):
        while True:
            print('Start parking and charging %d' % self.env.now)
            charge_duration = 5
            try:
                yield self.env.process(self.charge(charge_duration))
            except simpy.Interrupt:
                print('interrupted. Hope, that your battery is full enough...')
            
            print('Start driving at %d' % self.env.now)
            trip_duration = 2
            yield self.env.timeout(trip_duration)

    def charge(self, duration):
        yield self.env.timeout(duration)

def driver(env, car):
    yield env.timeout(3)
    car.action.interrupt()
    #yield env.process(driver(env,car)) ## when it shall interrupt every 3 timesteps

env = simpy.Environment()
car = Car(env)
env.process(driver(env,car))
#print(env.now) ## macht grundsätzlich den jetztingen Simulationszeitpunkt
env.run(until =30) #dieses run kommt aus der std_lib Simpy nicht von meiner methode


#
##Aleenas Beispiele
##-----------------------------------------------------------
#print("\n -->Aleenas Beispiele: \n")
#def main():
#	env = simpy.Environment() #tell process when environment they should be in
#	env.process(traffic_light(env)) # setup model, expecting to tell you the process that going to run in the environment, give the process  traffic_light
#	env.run(until = 120) # run the environment
#	print("simulation complete")
#
#
#def traffic_light(env): # the traffic light function needs to know that it is w/in this environment
#	while True:
#		print("Light turned GRN at t="+str(env.now)) #.now method gives you the simulation time (unitless)
#		yield env.timeout(30)
#		print("Light turned YELL at t="+str(env.now))
#		yield env.timeout(5)
#		print("Light turned RED at t="+str(env.now))
#		yield env.timeout(20)

#if __name__ == '__main__':
#	main()

