# PROJECT JULY-SEPTEMBRE 2019
# SOLVING THE N-BODIES PROBLEM / FUNCTIONS
# By Enguerran VIDAL

# This .py file contains the fucntions responsible for creating the different random distributions of massive objects
# we will later let evolve using our algorithms

###############################################################
#                           IMPORTS                           #
###############################################################

#-----------------MODULES
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib.animation

#-----------------PYTHON FILES
from __constants_conversions import*
from __functions import*
from __algorithms import*

###############################################################
#                           FONCTIONS                         #
###############################################################

def Random_Distribution_2D(distribution_type,distribution_seed,n_objects):
    ''' Creates a random distribution of objects in a 2D space by following a certain distribution type.
        Available 2D Distributions : Uniform Square, Uniform Disk, '''
    assert distribution_seed!=None and type(distribution_seed)==int
    np.random.seed(distribution_seed)
    n=n_objects
    objects=[]
    if distribution_type=='uniform_square':
        X=np.random.uniform(-10**(20), 10**(20),n)
        Y=np.random.uniform(-10**(20), 10**(20),n)
        dX=np.random.uniform(-10**(3), 10**(3),n)
        dY=np.random.uniform(-10**(3), 10**(3),n)
        M=np.random.uniform(10**(29),10**(32),n)
        for i in range(n):
            position=np.array([X[i],Y[i]])
            speed=np.array([dX[i],dY[i]])
            mass=M[i]
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_disk':
        for i in range(n):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            norm=np.sqrt(u*u+v*v)
            r=np.cbrt(np.random.random())
            U=np.array([u,v])
            position=(U/norm)*10**(20)*r # position on 1m sphere * radius 
            dx=np.random.uniform(-10**(3), 10**(3))
            dy=np.random.uniform(-10**(3), 10**(3))
            speed=np.array([dx,dy])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    return objects

def Random_Distribution_3D(distribution_type,distribution_seed,n_objects):
    ''' Creates a random distribution of objects in a 3D space by following a certain distribution type.
        Available 3D Distributions : Uniform Cube, Uniform Square, Uniform Sphere, Uniform Ball, Uniform Disk,
                                     Uniform Cylinder, Uniform Disk with cluster and Central Back Hole,
                                     Uniform Ball with cluster and Central Back Hole'''
    assert distribution_seed!=None and type(distribution_seed)==int
    np.random.seed(distribution_seed)
    n=n_objects
    objects=[]
    if distribution_type=='uniform_cube':
        X=np.random.uniform(-10**(20), 10**(20),n)
        Y=np.random.uniform(-10**(20), 10**(20),n)
        Z=np.random.uniform(-10**(20), 10**(20),n)
        dX=np.random.uniform(-10**(3), 10**(3),n)
        dY=np.random.uniform(-10**(3), 10**(3),n)
        dZ=np.random.uniform(-10**(3), 10**(3),n)
        M=np.random.uniform(10**(29),10**(32),n)
        for i in range(n):
            position=np.array([X[i],Y[i],Z[i]])
            speed=np.array([dX[i],dY[i],dZ[i]])
            mass=M[i]
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_square':
        X=np.random.uniform(-10**(20), 10**(20),n)
        Y=np.random.uniform(-10**(20), 10**(20),n)
        Z=np.zeros_like(X)
        dX=np.random.uniform(-10**(3), 10**(3),n)
        dY=np.random.uniform(-10**(3), 10**(3),n)
        dZ=np.zeros_like(dX)
        M=np.random.uniform(10**(29),10**(32),n)
        for i in range(n):
            position=np.array([X[i],Y[i],Z[i]])
            speed=np.array([dX[i],dY[i],dZ[i]])
            mass=M[i]
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_sphere':
        # We construct a sphere using the Muller method
        for i in range(n):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            w=np.random.normal(0,1)
            norm=np.sqrt(u*u+v*v+w*w)
            U=np.array([u,v,w])
            position=(U/norm)*10**(20) # position on 1m sphere * radius 
            dx=np.random.uniform(-10**(3), 10**(3))
            dy=np.random.uniform(-10**(3), 10**(3))
            dz=np.random.uniform(-10**(3), 10**(3))
            speed=np.array([dx,dy,dz])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_ball':
        # We construct a sphere using the Muller method
        for i in range(n):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            w=np.random.normal(0,1)
            norm=np.sqrt(u*u+v*v+w*w)
            r=np.cbrt(np.random.random())
            U=np.array([u,v,w])
            position=(U/norm)*10**(20)*r # position on 1m sphere * radius 
            dx=np.random.uniform(-10**(3), 10**(3))
            dy=np.random.uniform(-10**(3), 10**(3))
            dz=np.random.uniform(-10**(3), 10**(3))
            speed=np.array([dx,dy,dz])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_disk':
        for i in range(n):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            norm=np.sqrt(u*u+v*v)
            r=np.cbrt(np.random.random())
            U=np.array([u,v,0])
            position=(U/norm)*10**(20)*r # position on 1m disk * radius 
            dx=np.random.uniform(-10**(3), 10**(3))
            dy=np.random.uniform(-10**(3), 10**(3))
            dz=np.zeros_like(dx)
            speed=np.array([dx,dy,dz])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_cylinder':
        for i in range(n-1):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            z=np.random.uniform(-10**(20),10**(20))
            norm=np.sqrt(u*u+v*v)
            r=np.cbrt(np.random.random())
            U=np.array([u,v,0])
            position=(U/norm)*10**(20)*r # position on 1m disk * radius
            position[2]=z
            dx=np.random.uniform(-10**(3), 10**(3))
            dy=np.random.uniform(-10**(3), 10**(3))
            dz=np.zeros_like(dx)
            speed=np.array([dx,dy,dz])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
    if distribution_type=='uniform_disk_cluster+bh':
        for i in range(n-1):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            norm=np.sqrt(u*u+v*v)
            r=np.cbrt(np.random.random())
            U=np.array([u,v,0])
            position=(U/norm)*10**(20)*r # position on 1m disk * radius 
            dx=np.random.uniform(-10**(2), 10**(2))
            dy=np.random.uniform(-10**(2), 10**(2))
            dz=np.zeros_like(dx)
            speed=np.array([dx,dy,dz])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
        BH=Object(10**(36),np.array([0,0,0]),np.array([0,0,0]))
        objects.append(BH)
    if distribution_type=='uniform_ball_cluster+bh':
        # We construct a sphere using the Muller method
        for i in range(n-1):
            u=np.random.normal(0,1)
            v=np.random.normal(0,1)
            w=np.random.normal(0,1)
            norm=np.sqrt(u*u+v*v+w*w)
            r=np.cbrt(np.random.random())
            U=np.array([u,v,w])
            position=(U/norm)*10**(20)*r # position on 1m sphere * radius 
            dx=np.random.uniform(-10**(5), 10**(5))
            dy=np.random.uniform(-10**(5), 10**(5))
            dz=np.random.uniform(-10**(5), 10**(5))
            speed=np.array([dx,dy,dz])
            mass=np.random.uniform(10**(29),10**(32))
            objectI=Object(mass,position,speed)
            objects.append(objectI)
        BH=Object(10**(36),np.array([0,0,0]),np.array([0,0,0]))
        objects.append(BH)
        
    return objects

