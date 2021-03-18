# PROJECT JULY-SEPTEMBRE 2019
# SOLVING THE N-BODIES PROBLEM / FUNCTIONS
# By Enguerran VIDAL

# This file contains the multitude of functions used throughout this project, hence its importation in every single .py files


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


###############################################################
#                         FUNCTIONS                           #
###############################################################

def array_max_abs(array):
    ''' Works only for column vectors '''
    m=0
    for i in range(len(array)):
        x=abs(array[i])
        if x>m:
            m=x
    return m

def mean_value_array(X):
    ''' only works with column vectors of shape (n,)'''
    (n,)=X.shape
    S=0
    for i in range(n):
        S=S+X[i]
    X_bar=S/n
    return X_bar

def array_normal(X):
    return (X-min(X))/(max(X)-min(X))

#----------------------2D TRANSFORMATIONS

def points_vector_2D(pointI,pointJ):
    ''' Returns the distance between two points in 3D'''
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    return np.array([x,y])
    
def points_distance_2D(pointI,pointJ):
    ''' Returns the distance between two points in 3D'''
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    return np.sqrt(x**2+y**2)
        
def object_distance_2D(objectI,objectJ):
    ''' Returns the distance between two points in 3D'''
    v=object_vector_2D(objectI,objectJ)
    return np.sqrt(v[0]**2+v[1]**2)

def object_vector_2D(objectI,objectJ):
    ''' Returns the distance between two objects in 2D'''
    return objectJ.position-objectI.position

def vector_module_2D(vector):
    ''' Returns a vector's module in 2D'''
    return np.sqrt(vector[0]**2+vector[1]**2)

#----------------------3D TRANSFORMATIONS
def points_vector_3D(pointI,pointJ):
    ''' Returns the vector between two points in 3D'''
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    z=pointJ[2]-pointI[2]
    return np.array([x,y,z])
    
def points_distance_3D(pointI,pointJ):
    ''' Returns the distance between two points in 3D'''
    x=pointJ[0]-pointI[0]
    y=pointJ[1]-pointI[1]
    z=pointJ[2]-pointI[2]
    return np.sqrt(x**2+y**2+z**2)
        
def object_distance_3D(objectI,objectJ):
    ''' Returns the distance between two objects in 3D'''
    v=object_vector_3D(objectI,objectJ)
    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

def object_vector_3D(objectI,objectJ):
    ''' Returns the vector between two objects in 3D'''
    return objectJ.position-objectI.position

def vector_module_3D(vector):
    ''' Returns a vector's module in 3D'''
    return np.sqrt(vector[0]**2+vector[1]**2+vector[2]**2)


#------------------FILE READING

def translate_file(file_name):
    ''' Helps translate the formated parameters file in order to solve the NBody Problem ( is used in the main.py code )'''
    file=open(file_name,'r')
    lines=file.readlines()
    n=len(lines)
    parameters_labels=[None,None,None,None,None,None,None,None]
    parameters_values=[None,None,None,None,None,None,None,None]
    for i in range(n):
        lines[i]=lines[i].lower()
        content=lines[i].split()
        if len(content)>0:
            if content[0]=='dimension':
                parameters_labels[0]=content[0]
                parameters_values[0]=content[2]
            if content[0]=='algorithm' and content[1]=='method':
                parameters_labels[1]=content[0]+' '+content[1]
                parameters_values[1]=content[3]
            if content[0]=='distribution' and content[1]=='type':
                parameters_labels[2]=content[0]+' '+content[1]
                parameters_values[2]=content[3]
            if content[0]=='number' and content[1]=='of' and content[2]=='bodies':
                parameters_labels[3]=content[0]+' '+content[1]+' '+content[2]
                parameters_values[3]=int(content[4])
            if content[0]=='theta':
                parameters_labels[4]=content[0]
                parameters_values[4]=float(content[2])
            if content[0]=='frame' and content[1]=='length':
                parameters_labels[5]=content[0]+' '+content[1]
                parameters_values[5]=time_conversion(float(content[3]),content[4])
            if content[0]=='plot' and content[1]=='unit':
                parameters_labels[6]=content[0]+' '+content[1]
                parameters_values[6]=content[3]
            if content[0]=='distribution' and content[1]=='seed':
                parameters_labels[7]=content[0]+' '+content[1]
                parameters_values[7]=int(content[3])
    return parameters_labels,parameters_values


