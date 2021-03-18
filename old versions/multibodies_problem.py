#---------------IMPORTS

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib.animation

#---------------CLASSES

class Object():
    def __init__(self,mass,position,speed):
        self.mass=mass
        self.position=position
        self.speed=speed

class Simulation():
    def __init__(self,frame,Objects=None):
        self.frame=frame
        self.n_objects=0
        self.objects=[]
        if Objects!=None:
            n=len(Objects)
            for i in range(n):
                self.objects.append(Objects[i])
                self.n_objects=self.n_objects+1

    def add_object(self,Object):
        self.objects.append(Object)
        self.n_objects=self.n_objects+1
        
    def generate_objects(self,n,object_type=None):
        for i in range(n):
            position=np.array([random.uniform(-50, 50) for j in range(3)])
            speed=np.array([random.uniform(-0.05, 0.05) for j in range(3)])
            mass=random.uniform(10000000,10000000000)
            objectI=Object(mass,position,speed)
            self.objects.append(objectI)
            self.n_objects=self.n_objects+1
        
        
               
    def object_coords(self):
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        for i in range(self.n_objects):
            x=self.objects[i].position[0]
            y=self.objects[i].position[1]
            z=self.objects[i].position[2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
        return X,Y,Z
     
    def compute(self):
        new_positions=[]
        new_speeds=[]
        for i in range(self.n_objects):
            objectI=self.objects[i]
            acceleration=0
            for j in range(self.n_objects):
                if j!=i:
                    objectJ=self.objects[j]
                    m=objectJ.mass
                    d=object_distance(objectI,objectJ)
                    v=object_vector(objectI,objectJ)
                    a=(constant('G')*m*v)/(d**3)
                    acceleration=acceleration+a
            U=np.array([objectI.position,objectI.speed])
            New_U=U+self.frame*np.array([U[1],acceleration])
            new_positions.append(New_U[0])
            new_speeds.append(New_U[1])
        for i in range(self.n_objects):
            self.objects[i].position=new_positions[i]
            self.objects[i].speed=new_speeds[i]
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        for i in range(self.n_objects):
            x=new_positions[i][0]
            y=new_positions[i][1]
            z=new_positions[i][2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
        return X,Y,Z
         
    def run(self,duration=None):
        X,Y,Z=self.object_coords()
        print(X,Y,Z)
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        title = ax.set_title('3D Gravity') 
        graph = ax.scatter(X,Y,Z, c='darkblue', alpha=0.5)
        fig.show()
        continuer=1
        while continuer :
            plt.pause(0.0001)
            X,Y,Z=self.compute()
            graph._offsets3d = (X,Y,Z)
            plt.draw()


#---------------FONCTIONS

def object_distance(objectI,objectJ):
    v=object_vector(objectI,objectJ)
    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

def object_vector(objectI,objectJ):
    return objectJ.position-objectI.position

def constant(name):
    if name=='G':
        return 6.67430*10**(-11)


#---------------PROGRAM

S=Simulation(10)
S.generate_objects(2)
S.run()
