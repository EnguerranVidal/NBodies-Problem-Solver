# NBodies Problem Solver
 This project vows to achieve the creation of a physics simulation of the famous N-Bodies problem. The main Python class will choose a randomized distribution of massive bodies and will let them evolve, calculating their positions at the next time-step using the Barnes-Hut method in 2D or in 3D. This project's version is currently using Numpy elements but does not use it to calculate the next time-step positions. I will work on an optimized version in the future using Numpy and Cupy in order to get GPU-driven calculations ( because for now it is as slow as a snail...seriously). I do not recommend going over 500 simulated bodies at the same time as it quickly becomes laggy.

**Specifications** : This code uses Numpy, Matplotlib and Keyboard as its core libraries.


You will find the following files in this repository :

- **__algorithms.py** : It contains the Barnes-Hut and Standard Newtonian integration algorithms. Using the Euler method for next-step integration, it calculates the forces applied on the different massive bodies and calculate the next-step objects' positions.

- **__distributions.py** : It contains the functions responsible for creating a randomized distribution of massive bodies in a 2D or 3D space. A distribution type must be mentionned in order to get the function to use it.

- **__constants_conversions** : This .py file contains the conversion and constants calling functions used in the main.py file to aid vizualisation and help in displaying tangible time and space values/scales.

- **__functions.py** : This file contains the multitude of functions used in the rest of this project such as vector creations, distances and modules calculations and also the extraction of parameters from the parameters.txt file in the main.py project's main class.

- **main.py** : It contains the main N-Body Solver repsonsible for the plot display, keyboard event handling and the parameters being taken into account. This is the file to run in order to get the simulation running as it contains a small script.

- **parameters.txt** : Contains the main parameters directing the plot's appearance and speed as well as the Barnes-Hut algorithm accuracy and the number of massive bodies in our simulation.

I also put in a folder containing an even older version of this code made at the start of it, far before I learned to use GitHub.

Be aware that the plot can be controlled via mouse and the '+', '-', 'p' and 'q' keyboard keys. Once the simulation is upand running you can use :
- Your mouse as regular zoom and drag fucntions of the Matplotlib 2D and 3D plots.
- the '+' and '-' keyboard keys in order to increase or decrease the integration time-step value using a pseudo-logarithmic change of orders of magnitude.
- the 'q' key lets you end the simulation, I recommend pressing it for more than a second as the plt.draw() could reopen the simulation plot with nothing in it if not pressed for long enough.
- the 'p' key lets you pause the simulation if you want to better look around without the simulation still evolving.
