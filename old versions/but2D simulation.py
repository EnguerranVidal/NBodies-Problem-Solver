import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import matplotlib.animation


class Object():
    def __init__(self,mass,position,speed,index=None):
        self.mass=mass
        self.position=position
        self.speed=speed
        self.index=index


class QuadTree():
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
            d=points_distance(self.mass_center,objectI.position)
            if d!=0:
                t=width/d
                if t < theta or self.children[0]==None:
                    v=points_vector(objectI.position,self.mass_center)
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
                    
                 
            
class BHutSimulation2D():
                                    
    def __init__(self,frame,theta,Objects=None):
        self.frame=frame
        self.theta=theta
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
            position=np.array([random.uniform(-10**(20), 10**(20)) for j in range(2)])
            speed=np.array([random.uniform(-10**(6), 10**(6)) for j in range(2)])
            mass=random.uniform(10**(29),10**(40))
            objectI=Object(mass,position,speed)
            self.objects.append(objectI)
            self.n_objects=self.n_objects+1
            
    def object_coords(self):
        X=np.array([])
        Y=np.array([])
        for i in range(self.n_objects):
            x=self.objects[i].position[0]
            y=self.objects[i].position[1]
            X=np.append(X,x)
            Y=np.append(Y,y)
        return X,Y
    
    def norm_object_speed(self):
        dX=np.array([])
        dY=np.array([])
        for i in range(self.n_objects):
            dx=self.objects[i].position[0]
            dy=self.objects[i].position[1]
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
        M=np.sqrt(dX**2+dY**2)
        m=max(M)
        return M/m
    
    def object_speeds(self):
        dX=np.array([])
        dY=np.array([])
        for i in range(self.n_objects):
            dx=self.objects[i].position[0]
            dy=self.objects[i].position[1]
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
        return dX,dY   

    def compute(self):
        new_positions=[]
        new_speeds=[]
        X,Y=self.object_coords()
        m=max(array_max_abs(X),array_max_abs(Y))
        Q=QuadTree(self.objects,m,-m,m,-m)
        for i in range(self.n_objects):
            objectI=self.objects[i]
            acceleration=Q.forces(objectI,self.theta)
            U=np.array([objectI.position,objectI.speed])
            New_U=U+self.frame*np.array([U[1],acceleration])
            new_positions.append(New_U[0])
            new_speeds.append(New_U[1])
        for i in range(self.n_objects):
            self.objects[i].position=new_positions[i]
            self.objects[i].speed=new_speeds[i]
        X=[]
        Y=[]
        for i in range(self.n_objects):
            x=new_positions[i][0]
            y=new_positions[i][1]
            X.append(x)
            Y.append(y)
        return np.array([X,Y])
  

    def run(self):
        X,Y=self.object_coords()
        print(X,Y)
        M=self.norm_object_speed()
        print(M)
        fig = plt.figure()
        colormap='copper'
        ax = fig.add_subplot(111)
        title = ax.set_title('2D Gravity') 
        graph = ax.scatter(X,Y, c=M,cmap=colormap,vmin=0, vmax=1,edgecolor="k")
        plt.colorbar(graph, orientation="horizontal", pad=0.2)
        fig.show()
        continuer=1
        while continuer :
            plt.pause(0.0001)
            Offset=self.compute()
            Offset=Offset.T
            M=self.norm_object_speed()
            graph.set_offsets(Offset)
            graph.set_array(M)
            plt.draw()
        
        


def array_max_abs(array):
    ''' Works only for column vectors '''
    m=0
    for i in range(len(array)):
        x=abs(array[i])
        if x>m:
            m=x
    return m

def points_vector(pointI,pointJ):
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    return np.array([x,y])
    
def points_distance(pointI,pointJ):
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    return np.sqrt(x**2+y**2)
        
def object_distance(objectI,objectJ):
    v=object_vector(objectI,objectJ)
    return np.sqrt(v[0]**2+v[1]**2)

def object_vector(objectI,objectJ):
    return objectJ.position-objectI.position

def constant(name):
    if name=='G':
        return 6.67430*10**(-11)   

def years_to_seconds(t):
    return t*31536000

def hours_to_seconds(t):
    return t*3600

def minutes_to_seconds(t):
    return t*60

def days_to_seconds(t):
    return t*86400


S=BHutSimulation2D(years_to_seconds(5000),theta=1.25)
S.generate_objects(200)
S.run()
