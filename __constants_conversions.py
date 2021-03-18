# PROJECT JULY-SEPTEMBRE 2019
# SOLVING THE N-BODIES PROBLEM / FUNCTIONS
# By Enguerran VIDAL

# This .py file contains constants calling functions as well as regular conversions of time
# and distances in order to facilitate the mention and vizualisation of the results of our algorithms


###############################################################
#                         FUNCTIONS                           #
###############################################################

#---------------CONSTANTS

def constant(name):
    '''returns constnts from given names'''
    if name=='G':
        return 6.67430*10**(-11)


#---------------TIME CONVERSIONS

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

def time_stamp(t,stamp_type='barren'):
    time_value,time_unit=time(t)
    if stamp_type=='barren':
        return str(time_value)+' '+time_unit
    if stamp_type=='time_passed':
        return 'Time since t=0 : '+str(time_value)+' '+time_unit
    if stamp_type=='frame_length':
        return 'Frame Length : '+str(time_value)+' '+time_unit+'/f'

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


def reverse_time_conversion(t):
    '''returns t in a useful and meaningful unit
       t is initially in seconds'''
    if t>=Gy_to_seconds(1):
        return t/Gy_to_seconds(1),'Gy'
    else:
        if t>=My_to_seconds(1):
            return t/My_to_seconds(1),'My'
        else:
            if t>=Ky_to_seconds(1):
                return t/Ky_to_seconds(1),'Ky'
            else:
                if t>=years_to_seconds(1):
                    return t/years_to_seconds(1),'y'
                else:
                    if t>=days_to_seconds(1):
                        return t/days_to_seconds(1),'d'
                    else:
                        if t>=hours_to_seconds(1):
                            return t/hours_to_seconds(1),'h'
                        else:
                            if t>=minutes_to_seconds(1):
                                return t/minutes_to_seconds(1),'min'
                            else:
                                return t,'s'
                            
def time_conversion(value,unit):
    if unit=='Gy' or unit=='gy':
        return Gy_to_seconds(value)
    if unit=='My' or unit=='my':
        return My_to_seconds(value)
    if unit=='Ky' or unit=='ky':
        return Ky_to_seconds(value)
    if unit=='y':
        return years_to_seconds(value)
    if unit=='d':
        return days_to_seconds(value)
    if unit=='h':
        return hours_to_seconds(value)
    if unit=='min':
        return minutes_to_seconds(value)
    if unit=='s':
        return value

def unit_time_list():
    '''returns a list of usual time units'''
    return ['s','min','h','d','y','Ky','My','Gy']

def multiples(i):
    '''returns a useful multiple for scale units'''
    multiple_array=[[1,2,5,10,30],
                    [1,15,20,30],
                    [1,2,3,6,12],
                    [1,7,15,30,60,180],
                    [1,2,5,10,20,50,100,150,200,500],
                    [1,2,5,10,20,50,100,150,200,500],
                    [1,2,5,10,20,50,100,150,200,500],
                    [1,2,5,10,20,50,100,150,200,500]]
    return multiple_array[i]
    

#------------------DISTANCE CONVERSIONS

def ly_to_meters(length):
    '''transforms light years into meters'''
    return length*9.4607*10**(15)

def parsec_to_meters(length):
    '''transforms parsecs into meters'''
    return length*3.0857*10**(16)

def au_to_meters(length):
    '''transforms astronomical units into meters'''
    return length*149597870700

def Kparsec_to_meters(length):
    '''transforms Kpc into meters'''
    return length*parsec_to_meters(1000)

def Mparsec_to_meters(length):
    '''transforms Mpc units into meters'''
    return length*parsec_to_meters(1000000)

def Gparsec_to_meters(length):
    '''transforms Gpc into meters'''
    return length*parsec_to_meters(1000000000)

def distance_conversion(value,unit):
    if unit=='ly':
        return value/ly_to_meters(1)
    if unit=='pc':
        return value/parsec_to_meters(1)
    if unit=='au':
        return value/au_to_meters(1)
    if unit=='Kpc' or unit=='kpc':
        return value/Kparsec_to_meters(1)
    if unit=='m':
        return value
    if unit=='Mpc' or unit=='mpc':
        return value/Mparsec_to_meters(1)
    if unit=='Gpc' or unit=='gpc':
        return value/Gparsec_to_meters(1)
    if unit=='km':
        return value/1000

#------------------BASES CONVERSIONS

def spheric_to_cartesian(r,phi,theta):
    '''transforms spherical coordinates into cartesian cooridnates'''
    x=r*np.sin(theta)*np.cos(phi)
    y=r*np.sin(theta)*np.sin(phi)
    z=r*np.cos(theta)
    return x,y,z
