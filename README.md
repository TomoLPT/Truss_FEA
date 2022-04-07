# Truss_FEA

A python script to run a finite element analysis of a 3D pin jointed truss. Based on the stiffness method.

## Guide

Specify your parameters and geometry in the following format:
#### Nodes:
[(x, y, z, fixity_x, fixity_y, fixity_z), ...]  
where 1 = free, 0 = pinned for fixities  
example:  
(5, 1, 10, 0, 0, 1) ==> node at position x=5, y=1, z=10, pinned in the x and y direction but free to rotate in the z direction.

#### Forces:
[(node_index, x, y, z), ...]  
where x, y, z are magnitude of the forces in kN  
example:  
(3, 1, 1, 0) ==> force at node of index 3,  1kN in x-direction, 1kN in y-direction, 0 in z

#### Bars:
[(node_1, node_2), ...]  
example:  
(2, 5) ==> bar element from node at index 2 and 5  

#### Young's Modulus (E)
Only E value can be defined to be assigned to all elements  
Value in kN/m2

#### Bar area
Only one cross sectional area can be defined to be assigned to all elements  
Value in m2  

~~~python
nodes_test = [(0, 0, 0, 0, 0, 0), (5, 0, 0, 1, 1, 0), (0, 5, 0, 0, 0, 0)] #0 = restrained
forces_test = [(1, 0, -5, 0)] #node number, x_value, y_value
E_test = 200e6 #kN/m2
A_test = 100e-7
bars_test = [(0, 1), (1, 2)]
~~~
	
~~~sh
git clone https://github.com/TomoLPT/Truss_FEA
~~~
