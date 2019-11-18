![Schmetical Diagram](https://github.com/wangshaoyun/LJ_Fluid_Cell_list/blob/master/delete.jpg "Priciple Deletion")

# Monte Carlo Simulation of Lennard Jones Fluid by Cell Lists Method
## About the program
1. For the calculation of short-range potential, the most common strategy is the combination of Verlet list method and cell list method [1]. Or rather, use Verlet list to calcute force and use cell list method to update Verlet list.
2. In molecular dynamics, pure cell list method without Verlet list can directly used to calculte short-range force although it is litte slowly than the combination of Verlet list method and cell list method.
3. However, there is no corresponding cell list method in Monte Carlo. This program providea this program by doubly rotated linked lists technology [2,3].
4. This a program that simulate the equation of state and radial distribution function of LJ fluid by pure cell list method [2]. 
5. This program is veried through comparing the results with them by combination of Verlet list method and cell list method.
>[1] Frenkel D, Klein M, Parrrinello M. "Understanding Molecular Simulation: From Algorithm to Applications". 2nd ed. Singapore: Elsevier Pte Ltd, 2002.   
>[2] S. Y. Wang, C. H. Tong. Cell-lists Method for Monte Carlo Simulation, to be submitted.  
>[3] A. Drozdek, Data Structures and Algorithms in C++, Fourth ed., Cengage Learning, Boston, 2013.  
 
## About Parameter Files 
+ energy_data.txt: It is a parameter file which includes Lennard-Jones potential parameter and cutoff radius.
+ system_data.txt: It is a parameter file which includes system, time.  

## Compiling files
1. Determine the parameters such as box size, density, time step and so on.
2. Set the parameters in energy_data.txt and system_data.txt
3. Open terminal and enter the current file content
4. To compile the files:
```
$ make
```
  
5. To execute the files:
```
$ ./main &
```





