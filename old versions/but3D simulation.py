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


class OctTree():
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
            distance=points_distance(self.mass_center,objectI.position)
            t=width/distance
            if t < theta or self.children[0]==None:
                d=points_distance(self.mass_center,objectI.position)
                if d!=0:
                    v=points_vector(objectI.position,self.mass_center)
                    m=self.mass
                    a=(constant('G')*m*v)/(d**3)
                    return a
                else:
                    return np.array([0,0,0])
            else:
                a1=self.children[0].forces(objectI,theta)
                a2=self.children[1].forces(objectI,theta)
                a3=self.children[2].forces(objectI,theta)
                a4=self.children[3].forces(objectI,theta)
                a5=self.children[0].forces(objectI,theta)
                a6=self.children[1].forces(objectI,theta)
                a7=self.children[2].forces(objectI,theta)
                a8=self.children[3].forces(objectI,theta)
                return a1+a2+a3+a4+a5+a6+a7+a8
        else:
            return np.array([0,0,0])
                
            
class BHutSimulation3D():
                                    
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
        
    def generate_objects(self,n,distribution_type=None):
        if distribution_type=='cube':
            for i in range(n):
                position=np.array([random.uniform(-10**(16), 10**(16)) for j in range(3)])
                speed=np.array([random.uniform(-10**(3), 10**(3)) for j in range(3)])
                mass=random.uniform(10**(29),10**(32))
                objectI=Object(mass,position,speed)
                self.objects.append(objectI)
                self.n_objects=self.n_objects+1
        if distribution_type=='sphere':
            for i in range(n):
                radius=np.array([random.uniform(10**(-20), 10**(20))])
                phi=np.array([random.uniform(-np.pi, np.pi)])
                theta=np.array([random.uniform(0, np.pi)])
                x,y,z=spheric_to_cartesian(radius,phi,theta)
                position=np.array([x,y,z])
                position.shape=(3,)
                speed=np.array([random.uniform(-10**(3), 10**(3)) for j in range(3)])
                mass=random.uniform(10**(29),10**(33))
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

    def norm_object_speed(self):
        dX=np.array([])
        dY=np.array([])
        dZ=np.array([])
        for i in range(self.n_objects):
            dx=self.objects[i].position[0]
            dy=self.objects[i].position[1]
            dz=self.objects[i].position[2]
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
            dZ=np.append(dZ,dz)
        M=np.sqrt(dX**2+dY**2+dZ**2)
        return array_normal(M)
    
    def object_speeds(self):
        dX=np.array([])
        dY=np.array([])
        dZ=np.array([])
        for i in range(self.n_objects):
            dx=self.objects[i].position[0]
            dy=self.objects[i].position[1]
            dz=self.objects[i].position[2]
            dX=np.append(dX,dx)
            dY=np.append(dY,dy)
            dZ=np.append(dZ,dz)
        return dX,dY,dZ
        
    def compute(self):
        new_positions=[]
        new_speeds=[]
        X,Y,Z=self.object_coords()
        m=max([array_max_abs(X),array_max_abs(Y),array_max_abs(Z)])
        Q=OctTree(self.objects,m,-m,m,-m,m,-m)
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
  

    def run(self,speed_indicator=True,unit='m'):
        if speed_indicator==True:
            X,Y,Z=self.object_coords()
            if unit=='ly':
                X,Y,Z=X/ly_to_meters(1),Y/ly_to_meters(1),Z/ly_to_meters(1)
            t=0
            M=self.norm_object_speed()
            colormap='copper'
            fig = plt.figure()
            ax = fig.gca(projection='3d')
            title = ax.set_title('3D Gravity') 
            graph = ax.scatter(X,Y,Z, c=M,cmap=colormap,vmin=0, vmax=1,edgecolor="k")
            plt.colorbar(graph, orientation="horizontal", pad=0.2)
            props = dict(boxstyle='round',facecolor='wheat',alpha=0.5)
            text_content=time_stamp(t,'time-passed')
            time_text=fig.text(0.05,-0.2,text_content,transform=ax.transAxes,fontsize=10,verticalalignment='top',bbox=props)
            fig.show()
            continuer=1
            while continuer :
                plt.pause(0.0001)
                X,Y,Z=self.compute()
                if unit=='ly':
                    X,Y,Z=X/ly_to_meters(1),Y/ly_to_meters(1),Z/ly_to_meters(1)
                M=self.norm_object_speed()
                graph._offsets3d = (X,Y,Z)
                graph.set_array(M)
                t=t+self.frame
                time_text.set_text(time_stamp(t,'time-passed'))
                plt.draw()
        else:
            X,Y,Z=self.object_coords()
            if unit=='ly':
                X,Y,Z=X/ly_to_meters(1),Y/ly_to_meters(1),Z/ly_to_meters(1)
            t=0
            colormap='copper'
            fig = plt.figure()
            colors=np.full(X.shape,1)
            ax = fig.gca(projection='3d')
            title = ax.set_title('3D Gravity') 
            graph = ax.scatter(X,Y,Z, c=colors,cmap=colormap,vmin=0, vmax=1,edgecolor="k")
            props = dict(boxstyle='round',facecolor='wheat',alpha=0.5)
            text_content=time_stamp(t,'time-passed')
            time_text=fig.text(0.05,-0.2,text_content,transform=ax.transAxes,fontsize=10,verticalalignment='top',bbox=props)
            fig.show()
            continuer=1
            while continuer :
                plt.pause(0.0001)
                X,Y,Z=self.compute()
                if unit=='ly':
                    X,Y,Z=X/ly_to_meters(1),Y/ly_to_meters(1),Z/ly_to_meters(1)
                graph._offsets3d = (X,Y,Z)
                t=t+self.frame
                time_text.set_text(time_stamp(t,'time-passed'))
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
    z=pointJ[2]-pointI[2]
    return np.array([x,y,z])
    
def points_distance(pointI,pointJ):
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    z=pointJ[2]-pointI[2]
    return np.sqrt(x**2+y**2+z**2)
        
def object_distance(objectI,objectJ):
    v=object_vector(objectI,objectJ)
    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

def object_vector(objectI,objectJ):
    return objectJ.position-objectI.position

def vector_module(vector):
    return np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)

def array_normal(X):
    return (X-min(X))/(max(X)-min(X))
             
def mean_value_array(X):
    ''' only works with column vectors of shape (n,)'''
    (n,)=X.shape
    S=0
    for i in range(n):
        S=S+X[i]
    X_bar=S/n
    return X_bar

def constant(name):
    if name=='G':
        return 6.67430*10**(-11)

def spheric_to_cartesian(r,phi,theta):
    x=r*np.sin(theta)*np.cos(phi)
    y=r*np.sin(theta)*np.sin(phi)
    z=r*np.cos(theta)
    return x,y,z
    

def Gy_to_seconds(t):
    return years_to_seconds(t)*1000000000
    
def My_to_seconds(t):
    return years_to_seconds(t)*1000000

def Ky_to_seconds(t):
    return years_to_seconds(t)*1000

def years_to_seconds(t):
    return t*31536000

def hours_to_seconds(t):
    return t*3600

def minutes_to_seconds(t):
    return t*60

def days_to_seconds(t):
    return t*86400

def ly_to_meters(length):
    return length*9.4607*10**(15)

def time_stamp(t,stamp_type='barren'):
    time_value,time_unit=time(t)
    if stamp_type=='barren':
        return str(time_value)+' '+time_unit
    if stamp_type=='time-passed':
        return 'Time since t=0 : '+str(time_value)+' '+time_unit

def time(t):
    if t>=Gy_to_seconds(1):
        return round(t/Gy_to_seconds(1),2),'G years'
    else:
        if t>=My_to_seconds(1):
            return round(t/My_to_seconds(1),2),'M years'
        else:
            if t>=Ky_to_seconds(1):
                return round(t/Ky_to_seconds(1),2),'K years'
            else:
                if t>=years_to_seconds(1):
                    return round(t/years_to_seconds(1),2),'years'
                else:
                    if t>=days_to_seconds(1):
                        return round(t/days_to_seconds(1),2),'days'
                    else:
                        if t>=hours_to_seconds(1):
                            return round(t/hours_to_seconds(1),2),'hours'
                        else:
                            if t>=minutes_to_seconds(1):
                                return round(t/minutes_to_seconds(1),2),'minutes'
                            else:
                                return round(t,2),'seconds'

S=BHutSimulation3D(Ky_to_seconds(100),theta=1.25)
S.generate_objects(2000,distribution_type='sphere')
S.run(speed_indicator=False,unit='ly')
