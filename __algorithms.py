# PROJECT JULY-SEPTEMBRE 2019
# SOLVING THE N-BODIES PROBLEM / CLASSES
# By Enguerran VIDAL

# This file contains the different kind of engines used further in main.py to
# calculate the state of our NBody problem at the next time-step.

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

###############################################################
#                           CLASSES                           #
###############################################################


#------------OBJECTS
class Object():
    ''' A class creating a point in space with mass '''
    def __init__(self,mass,position,speed,index=None):
        self.mass=mass
        self.position=position
        self.speed=speed
        self.index=index


#------------SUB ALGORITHMS
class OctTree():
    ''' A class creating an object permiting the division of 3D space in subdividing cubes
to finally use it to compute a Barnes-Hut algorithm ( solving the N-bodies problem ).'''
    def __init__(self,objects,a,b,c,d,e,f):
        self.objects=objects
        self.n_objects=len(self.objects)
        self.children=[None,None,None,None,None,None,None,None]
        self.a=a
        self.b=b
        self.c=c
        self.d=d
        self.e=e
        self.f=f
        self.m_center()
        self.build()

    def build(self):
        if self.n_objects>1:
            TNW=[]
            TNE=[]
            TSW=[]
            TSE=[]
            BNW=[]
            BNE=[]
            BSW=[]
            BSE=[]
            for i in range(self.n_objects):
                m1=(self.a+self.b)/2
                m2=(self.c+self.d)/2
                m3=(self.e+self.f)/2
                if self.objects[i].position[0]>m1:
                    if self.objects[i].position[1]>m2:
                        if self.objects[i].position[2]>m3:
                            TNE.append(self.objects[i])
                        else:
                            BNE.append(self.objects[i])
                    else :
                        if self.objects[i].position[2]>m3:
                            TSE.append(self.objects[i])
                        else:
                            BSE.append(self.objects[i])
                else :
                    if self.objects[i].position[1]>m2:
                        if self.objects[i].position[2]>m3:
                            TNW.append(self.objects[i])
                        else:
                            BNW.append(self.objects[i])
                    else :
                        if self.objects[i].position[2]>m3:
                            TSW.append(self.objects[i])
                        else:
                            BSW.append(self.objects[i])
                        
            self.children[0]=OctTree(TNW,m1,self.b,self.c,m2,self.e,m3)
            self.children[1]=OctTree(TNE,self.a,m1,self.c,m2,self.e,m3)
            self.children[2]=OctTree(TSW,m1,self.b,m2,self.d,self.e,m3)
            self.children[3]=OctTree(TSE,self.a,m1,m2,self.d,self.e,m3)
            self.children[4]=OctTree(BNW,m1,self.b,self.c,m2,m3,self.f)
            self.children[5]=OctTree(BNE,self.a,m1,self.c,m2,m3,self.f)
            self.children[6]=OctTree(BSW,m1,self.b,m2,self.d,m3,self.f)
            self.children[7]=OctTree(BSE,self.a,m1,m2,self.d,m3,self.f)

    def m_center(self):
        if self.n_objects>0:
            xc=0
            yc=0
            zc=0
            m=0
            for i in range(self.n_objects):
                xc=xc+self.objects[i].mass*self.objects[i].position[0]
                yc=yc+self.objects[i].mass*self.objects[i].position[1]
                yc=yc+self.objects[i].mass*self.objects[i].position[2]
                m=m+self.objects[i].mass
            self.mass_center=[xc/m,yc/m,zc/m]
            self.mass=m
            
    def forces(self,objectI,theta):
        if self.n_objects>0:
            width=self.a-self.b
            d=points_distance_3D(self.mass_center,objectI.position)
            if d!=0:
                t=width/d
                if t < theta or self.children[0]==None:
                    v=points_vector_3D(objectI.position,self.mass_center)
                    m=self.mass
                    a=(constant('G')*m*v)/(d**3)
                    return a
                else:
                    a1=self.children[0].forces(objectI,theta)
                    a2=self.children[1].forces(objectI,theta)
                    a3=self.children[2].forces(objectI,theta)
                    a4=self.children[3].forces(objectI,theta)
                    a5=self.children[4].forces(objectI,theta)
                    a6=self.children[5].forces(objectI,theta)
                    a7=self.children[6].forces(objectI,theta)
                    a8=self.children[7].forces(objectI,theta)
                    return a1+a2+a3+a4
            else:
                return np.array([0,0,0])
        else:
            return np.array([0,0,0])


class QuadTree():
    ''' A class creating an object permiting the division of 2D space in subdividing squares
to finally use it to compute a Barnes-Hut algorithm ( solving the N-bodies problem ).'''
    def __init__(self,objects,a,b,c,d):
        self.objects=objects
        self.n_objects=len(self.objects)
        self.children=[None,None,None,None]
        self.a=a
        self.b=b
        self.c=c
        self.d=d
        self.m_center()
        self.build()

    def build(self):
        if self.n_objects>1:
            NW=[]
            NE=[]
            SW=[]
            SE=[]
            for i in range(self.n_objects):
                m1=(self.a+self.b)/2
                m2=(self.c+self.d)/2
                if self.objects[i].position[0]>m1:
                    if self.objects[i].position[1]>m2:
                        NE.append(self.objects[i])
                    else :
                        SE.append(self.objects[i])
                else :
                    if self.objects[i].position[1]>m2:
                        NW.append(self.objects[i])
                    else :
                        SW.append(self.objects[i])
            self.children[0]=QuadTree(NW,m1,self.b,self.c,m2)
            self.children[1]=QuadTree(NE,self.a,m1,self.c,m2)
            self.children[2]=QuadTree(SW,m1,self.b,m2,self.d)
            self.children[3]=QuadTree(SE,self.a,m1,m2,self.d)

    def m_center(self):
        if self.n_objects>0:
            xc=0
            yc=0
            m=0
            for i in range(self.n_objects):
                xc=xc+self.objects[i].mass*self.objects[i].position[0]
                yc=yc+self.objects[i].mass*self.objects[i].position[1]
                m=m+self.objects[i].mass
            self.mass_center=[xc/m,yc/m]
            self.mass=m
                   
    def forces(self,objectI,theta):
        if self.n_objects>0:
            width=self.a-self.b
            d=points_distance_2D(self.mass_center,objectI.position)
            if d!=0:
                t=width/d
                if t < theta or self.children[0]==None:
                    v=points_vector_2D(objectI.position,self.mass_center)
                    m=self.mass
                    a=(constant('G')*m*v)/(d**3)
                    return a
                else:
                    a1=self.children[0].forces(objectI,theta)
                    a2=self.children[1].forces(objectI,theta)
                    a3=self.children[2].forces(objectI,theta)
                    a4=self.children[3].forces(objectI,theta)
                    return a1+a2+a3+a4
            else:
                return np.array([0,0])
        else:
            return np.array([0,0])


#----------------SIMULATION_ALGORITHMS

class BHut_Algorithm_3D():
    ''' Class using the Octree class above in order to solve the N-Bodies problem in 3D using the Barnes-Hut Algorithm'''                                 
    def __init__(self,frame,theta):
        self.frame=frame
        self.theta=theta

    def objects_positions(self,Objects):
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        n=len(Objects)
        for i in range(n):
            x=Objects[i].position[0]
            y=Objects[i].position[1]
            z=Objects[i].position[2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
        return X,Y,Z
    
    def update_frame(new_frame):
        self.frame=new_frame
        
    def compute(self,Objects):
        new_positions=[]
        new_speeds=[]
        n=len(Objects)
        X,Y,Z=self.objects_positions(Objects)
        m=max([array_max_abs(X),array_max_abs(Y),array_max_abs(Z)])
        Q=OctTree(Objects,m,-m,m,-m,m,-m)
        for i in range(n):
            objectI=Objects[i]
            acceleration=Q.forces(objectI,self.theta)
            U=np.array([objectI.position,objectI.speed])
            New_U=U+self.frame*np.array([U[1],acceleration])
            new_positions.append(New_U[0])
            new_speeds.append(New_U[1])
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        dX=np.array([])
        dY=np.array([])
        dZ=np.array([])
        for i in range(n):
            x=new_positions[i][0]
            y=new_positions[i][1]
            z=new_positions[i][2]
            dx=new_speeds[i][0]
            dy=new_speeds[i][1]
            dz=new_speeds[i][2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
            dZ=np.append(dZ,dz)
        return X,Y,Z,dX,dY,dZ


class BHut_Algorithm_2D():
    ''' Class using the Quadtree class above in order to solve the N-Bodies problem in 2D using the Barnes-Hut Algorithm'''                                 
    def __init__(self,frame,theta):
        self.frame=frame
        self.theta=theta

    def objects_positions(self,Objects):
        X=np.array([])
        Y=np.array([])
        n=len(Objects)
        for i in range(n):
            x=Objects[i].position[0]
            y=Objects[i].position[1]
            X=np.append(X,x)
            Y=np.append(Y,y)
        return X,Y

    def update_frame(new_frame):
        self.frame=new_frame

    def compute(self,Objects):
        new_positions=[]
        new_speeds=[]
        X,Y=self.objects_positions(Objects)
        m=max(array_max_abs(X),array_max_abs(Y))
        Q=QuadTree(Objects,m,-m,m,-m)
        n=len(Objects)
        for i in range(n):
            objectI=Objects[i]
            acceleration=Q.forces(objectI,self.theta)
            U=np.array([objectI.position,objectI.speed])
            New_U=U+self.frame*np.array([U[1],acceleration])
            new_positions.append(New_U[0])
            new_speeds.append(New_U[1])
        X=np.array([])
        Y=np.array([])
        dX=np.array([])
        dY=np.array([])
        for i in range(n):
            x=new_positions[i][0]
            y=new_positions[i][1]
            dx=new_speeds[i][0]
            dy=new_speeds[i][1]
            X=np.append(X,x)
            Y=np.append(Y,y)
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
        return X,Y,dX,dY

class Standard_Algorithm_3D():
    ''' Class solving the N-Bodies problem in 3D using a regular integrator algorithm without any real approximations.''' 
    def __init__(self,frame):
        self.frame=frame    
               
    def objects_positions(self,Objects):
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        n=len(Objects)
        for i in range(n):
            x=Objects[i].position[0]
            y=Objects[i].position[1]
            z=Objects[i].position[2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
        return X,Y,Z

    def update_frame(new_frame):
        self.frame=new_frame
     
    def compute(self,Objects):
        new_positions=[]
        new_speeds=[]
        n=len(Objects)
        for i in range(n):
            objectI=Objects[i]
            acceleration=0
            for j in range(n):
                if j!=i:
                    objectJ=Objects[j]
                    m=objectJ.mass
                    d=object_distance_3D(objectI,objectJ)
                    v=object_vector_3D(objectI,objectJ)
                    a=(constant('G')*m*v)/(d**3)
                    acceleration=acceleration+a
            U=np.array([objectI.position,objectI.speed])
            New_U=U+self.frame*np.array([U[1],acceleration])
            new_positions.append(New_U[0])
            new_speeds.append(New_U[1])
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        dX=np.array([])
        dY=np.array([])
        dZ=np.array([])
        for i in range(n):
            x=new_positions[i][0]
            y=new_positions[i][1]
            z=new_positions[i][2]
            dx=new_speeds[i][0]
            dy=new_speeds[i][1]
            dz=new_speeds[i][2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
            dZ=np.append(dZ,dz)
        return X,Y,Z,dX,dY,dZ


class Standard_Algorithm_2D():
    '''Class solving the N-Bodies problem in 2D using a regular integrator algorithm without any real approximations.'''
    def __init__(self,frame):
        self.frame=frame    
               
    def objects_positions(self,Objects):
        X=np.array([])
        Y=np.array([])
        n=len(Objects)
        for i in range(n):
            x=Objects[i].position[0]
            y=Objects[i].position[1]
            X=np.append(X,x)
            Y=np.append(Y,y)
        return X,Y
    
    def update_frame(new_frame):
        self.frame=new_frame
     
    def compute(self,Objects):
        new_positions=[]
        new_speeds=[]
        n=len(Objects)
        for i in range(n):
            objectI=Objects[i]
            acceleration=0
            for j in range(n):
                if j!=i:
                    objectJ=Objects[j]
                    m=objectJ.mass
                    d=object_distance_2D(objectI,objectJ)
                    v=object_vector_2D(objectI,objectJ)
                    a=(constant('G')*m*v)/(d**3)
                    acceleration=acceleration+a
            U=np.array([objectI.position,objectI.speed])
            New_U=U+self.frame*np.array([U[1],acceleration])
            new_positions.append(New_U[0])
            new_speeds.append(New_U[1])
        X=[]
        Y=[]
        dX=[]
        dY=[]
        for i in range(self.n_objects):
            x=new_positions[i][0]
            y=new_positions[i][1]
            dx=new_speeds[i][0]
            dy=new_speeds[i][1]
            X.append(x)
            Y.append(y)
            dX.append(dx)
            dY.append(dy)
        return X,Y,dX,dY
