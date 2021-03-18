# PROJECT JULY-SEPTEMBRE 2019
# SOLVING THE N-BODIES PROBLEM / FUNCTIONS
# By Enguerran VIDAL

# This file contains the main class of this project.

###############################################################
#                           IMPORTS                           #
###############################################################

#-----------------MODULES
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import keyboard

#-----------------PYTHON FILES
from __constants_conversions import*
from __algorithms import*
from __functions import*
from __distributions import*

###############################################################
#                        MAIN CLASS                           #
###############################################################

class N_Bodies_Problem():
    ''' Main class, taking advantage of the other project's classes. Its goal is to help the user choose its parameters as well as import them from a text file formated in a certain way.
It will then plot the objects position and update it as the algorithm calculates their places at the next time step. It also gives the user the ability to speed up the simulation, for the price
of further error propagation in the simulation itself, meaning the objects' positions will stray further and further from the truth.'''
    def __init__(self,parameters_file):
        self.Objects=[]
        # Distributing the parameters from the file given after tra
        parameters_labels,parameters_values=translate_file(parameters_file)
        print(parameters_labels)
        print(parameters_values)
        self.dimension=parameters_values[0]
        self.algorithm_method=parameters_values[1]
        self.distribution_type=parameters_values[2]
        self.number_of_bodies=parameters_values[3]
        self.theta=parameters_values[4]
        self.frame=parameters_values[5]
        self.original_frame=parameters_values[5]
        self.plot_unit=parameters_values[6]
        self.distribution_seed=parameters_values[7]
        # Redirecting infos to the different choices given
        self.distribution_choosing()
        self.method_choosing()
        
    def method_choosing(self):
        if self.dimension=='2d':
            if self.algorithm_method=='bhut':
                self.Algorithm=BHut_Algorithm_2D(self.frame,self.theta)
            if self.algorithm_method=='standard':
                self.Algorithm=Standard_Algorithm_2D(self.frame)
        if self.dimension=='3d':
            if self.algorithm_method=='bhut':
                self.Algorithm=BHut_Algorithm_3D(self.frame,self.theta)
            if self.algorithm_method=='standard':
                self.Algorithm=Standard_Algorithm_3D(self.frame)
        print(self.Algorithm)

    def distribution_choosing(self):
        if self.dimension=='2d':
            self.Objects=Random_Distribution_2D(self.distribution_type,self.distribution_seed,self.number_of_bodies)
        if self.dimension=='3d':
            self.Objects=Random_Distribution_3D(self.distribution_type,self.distribution_seed,self.number_of_bodies)

    def objects_update_3D(self,X,Y,Z,dX,dY,dZ):
        n=len(self.Objects)
        for i in range(n):
            self.Objects[i].position=np.array([X[i],Y[i],Z[i]])
            self.Objects[i].speed=np.array([dX[i],dY[i],dZ[i]])
            
    def objects_update_2D(self,X,Y,dX,dY):
        n=len(self.Objects)
        for i in range(n):
            self.Objects[i].position=np.array([X[i],Y[i]])
            self.Objects[i].speed=np.array([dX[i],dY[i]])
            
    def objects_coordinates_2D(self):
        X=np.array([])
        Y=np.array([])
        n=len(self.Objects)
        for i in range(n):
            x=self.Objects[i].position[0]
            y=self.Objects[i].position[1]
            X=np.append(X,x)
            Y=np.append(Y,y)
        return X,Y
    
    def objects_coordinates_3D(self):
        X=np.array([])
        Y=np.array([])
        Z=np.array([])
        n=len(self.Objects)
        for i in range(n):
            x=self.Objects[i].position[0]
            y=self.Objects[i].position[1]
            z=self.Objects[i].position[2]
            X=np.append(X,x)
            Y=np.append(Y,y)
            Z=np.append(Z,z)
        return X,Y,Z

    def frame_updating(self,sign):
        if sign=='+': # Goes to a higher multiple of time
            time_value,time_unit=reverse_time_conversion(self.frame)
            unit_list= unit_time_list()
            i=unit_list.index(time_unit)
            coefficients=multiples(i)
            if time_value in coefficients:
                j=coefficients.index(time_value)
                if j==len(coefficients)-1:
                    if i==len(unit_list)-1:
                        print("WARNING : Reached maximum simulation speed ! Can't increase frame length further.")
                    else :
                        time_value=multiples(i+1)[0]
                        time_unit=unit_list[i+1]
                        self.frame=time_conversion(time_value,time_unit)
                else :
                    time_value=coefficients[j+1]
                    self.frame=time_conversion(time_value,time_unit)
            else:
                if time_value>coefficients[-1]:
                    if i==len(unit_list)-1:
                        print("WARNING : Reached maximum simulation speed ! Can't increase frame length further.")
                    else :
                        time_value=multiples(i+1)[0]
                        time_unit=unit_list[i+1]
                        self.frame=time_conversion(time_value,time_unit)
                else:
                    j=0
                    while time_value<coefficients[j]:
                        j=j+1
                    time_value=coefficients[j]
                    self.frame=time_conversion(time_value,time_unit)
                    
        if sign=='-': # Goes to a lower multiple of time
            time_value,time_unit=reverse_time_conversion(self.frame)
            unit_list= unit_time_list()
            i=unit_list.index(time_unit)
            coefficients=multiples(i)
            if time_value in coefficients:
                j=coefficients.index(time_value)
                if j==0:
                    if i==0:
                        print("WARNING : Reached minimum simulation speed ! Can't decrease frame length further.")
                    else :
                        time_value=multiples(i-1)[-1]
                        time_unit=unit_list[i-1]
                        self.frame=time_conversion(time_value,time_unit)
                else :
                    time_value=coefficients[j-1]
                    self.frame=time_conversion(time_value,time_unit)
            else :
                if time_value<coefficients[0]:
                    if i==0:
                        print("WARNING : Reached minimum simulation speed ! Can't decrease frame length further.")
                    else :
                        time_value=multiples(i-1)[-1]
                        time_unit=unit_list[i-1]
                        self.frame=time_conversion(time_value,time_unit)
                else :
                    j=len(coefficients)-1
                    while time_value>coefficients[j]:
                        j=j-1
                    time_value=coefficients[j]
                    self.frame=time_conversion(time_value,time_unit)

        if sign=='*': # Change the frame length to its original value
            self.frame=self.original_frame


    def run(self):
        if self.dimension=='2d':
            t=0
            # Defining plot
            fig = plt.figure()
            ax = fig.add_subplot(111)
            title = ax.set_title('2D Gravity')
            X,Y=self.objects_coordinates_2D()
            X=distance_conversion(X,self.plot_unit)
            Y=distance_conversion(Y,self.plot_unit)
            graph = ax.scatter(X,Y, c='y',edgecolor="k")
            # Defining the time label
            props = dict(boxstyle='round',facecolor='wheat',alpha=0.5)
            text_content=time_stamp(t,'time-passed')
            time_text=fig.text(0.05,0.05,text_content,transform=ax.transAxes,fontsize=10,verticalalignment='top',bbox=props)
            # Plotting everything
            fig.show()
            # Main loop
            continuer=1
            while continuer :
                plt.pause(0.01)
                X,Y,dX,dY=self.Algorithm.compute(self.Objects)
                self.objects_update_2D(X,Y,dX,dY)
                # Conversion of the coordinates to plot_unit
                X=distance_conversion(X,self.plot_unit)
                Y=distance_conversion(Y,self.plot_unit)
                # Updating the scatter plot
                Offset=np.array([X,Y])
                Offset=Offset.T
                graph.set_offsets(Offset)
                # Updting the time label
                t=t+self.frame
                time_text.set_text(time_stamp(t,'time-passed'))
                # Drawing the changes
                plt.draw()
        if self.dimension=='3d':
            t=0
            print(len(self.Objects))
            # Defining plot
            fig = plt.figure(figsize=(6,6))
            ax1 = fig.gca(projection='3d')
            title = ax1.set_title('3D Gravity')
            X,Y,Z=self.objects_coordinates_3D()
            # Conversion of the coordinates to plot_unit
            X=distance_conversion(X,self.plot_unit)
            Y=distance_conversion(Y,self.plot_unit)
            Z=distance_conversion(Z,self.plot_unit)
            graph = ax1.scatter(X,Y,Z, c='y',edgecolor="k")
            # Defining the time and frame label
            props = dict(boxstyle='round',facecolor='wheat',alpha=0.5)
            text_content=time_stamp(t,'time-passed')
            time_text1=fig.text(0.05,0.05,text_content,transform=ax1.transAxes,fontsize=10,verticalalignment='top',bbox=props)
            text_content2=time_stamp(self.frame,'frame_length')
            time_text2=fig.text(0.05,0.05,text_content2,transform=ax1.transAxes,fontsize=10,verticalalignment='top',bbox=props)
            # Plotting everything
            fig.show()
            plt.grid()
            # Main loop
            continue_factor=1
            while continue_factor :
                if keyboard.is_pressed('spacebar'):  # if key 'space' is pressed 
                    print('Pausing!')
                    plt.pause(1)
                    while True:
                        plt.pause(0.0001)
                        if keyboard.is_pressed('spacebar'):  # if key 'space' is pressed 
                            print('Depausing!')
                            break
                        if keyboard.is_pressed('q'): # if key 'q' is pressed 
                            continue_factor=0
                            break
                        if keyboard.is_pressed('+'):  # if key '+' is pressed
                            plt.pause(0.001)
                            self.frame_updating('+')
                        if keyboard.is_pressed('-'):# if key '-' is pressed 
                            plt.pause(0.001)
                            self.frame_updating('-')
                        if keyboard.is_pressed('*'):  # if key '*' is pressed 
                            plt.pause(0.001)
                            self.frame_updating('*')
                if keyboard.is_pressed('+'):  # if key '+' is pressed 
                    plt.pause(0.001)
                    self.frame_updating('+')
                if keyboard.is_pressed('-'):  # if key '-' is pressed 
                    plt.pause(0.001)
                    self.frame_updating('-')
                if keyboard.is_pressed('*'):  # if key '*' is pressed 
                    plt.pause(0.001)
                    self.frame_updating('*')
                if keyboard.is_pressed('q') or continue_factor==0:  # if key 'q' is pressed 
                    print('Quitting!')
                    break
                plt.pause(0.001)
                self.Algorithm.frame=self.frame
                X,Y,Z,dX,dY,dZ=self.Algorithm.compute(self.Objects)
                self.objects_update_3D(X,Y,Z,dX,dY,dZ)
                # Conversion of the coordinates to plot_unit
                X=distance_conversion(X,self.plot_unit)
                Y=distance_conversion(Y,self.plot_unit)
                Z=distance_conversion(Z,self.plot_unit)
                # Updating the scatter plot
                graph._offsets3d = (X,Y,Z)
                # Updting the time label
                t=t+self.frame
                time_text1.set_text(time_stamp(t,'time-passed'))
                time_text2.set_text(time_stamp(self.frame,'frame_length'))
                # Drawing the changes
                plt.draw()
            plt.close()


###############################################################
#                       MAIN PROGRAM                          #
###############################################################

Simulation=N_Bodies_Problem('parameters.txt')
Simulation.run()
