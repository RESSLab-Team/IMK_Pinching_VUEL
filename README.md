---------------
DESCRIPTION																						   
---------------

VUEL SUBROUTINE FOR IMK SPRING WITH PINCHING/PEAK-ORIENTED BEHAVIOUR									
																										
This file provides a guide for the implementation of the user-defined element (VUEL) VU1.for in an Abaqus input
file. The element may be used in Abaqus/Explicit only.

The element consists of two uncoupled springs in two directions: a non-linear spring and an elastic spring.The 
non-linear spring simulates the Ibarra-Medina-Krawinkler deterioration model with peak-oriented or pinching
hysteretic response. The behavior incorporates basic strength, post-capping strength, accelerated reloading 
stiffness, and unloading stiffness deterioration. Details of the model are found in Ibarra et al. (2005).


------------------
IMPLEMENTATION																						   
------------------

Below is an example of the implementation in an Abaqus input file. The element corresponds to a zero-length
spring that represents the behavior of a shear stud (El Jisr et al. 2020).

The properties of the element are as follows:


++Monotonic Backbone Properties++
- Elastic stiffness, Ke = 92000
- Pre-capping plastic deformation in the positive loading direction, Up_pos = 6.6
- Post-capping plastic deformation in the positive loading direction, Upc_pos = 11.0
- Effective yield force the positive loading direction, Fy_pos = 76000
- Maximum-to-effective yield force ratio in the positive loading direction, Fmax_Fy_pos = 1.07
- Residual-to-effective yield force ratio in the  positive loading direction, res_pos = 0.2
- Ultimate deformation in the positive loading direction, Uu_pos = 15.0
- Pre-capping plastic deformation in the negative loading direction, Up_neg = 10.2 
- Post-capping plastic deformation in the negative loading direction, Upc_neg = 5.0
- Effective yield force the negative loading direction, Fy_neg = 34000
- Maximum-to-effective yield force ratio in the positive loading direction, Fmax_Fy_neg = 1.07
- Residual-to-effective yield force ratio in the  positive loading direction, res_neg = 0.2
- Ultimate deformation in the positive loading direction, Uu_neg = 15.0


++Cyclic Deterioration Properties++
- Cyclic deterioration parameter for strength deterioration, lambda_S = 40.0
- Cyclic deterioration parameter for post-capping strength deterioration, lambda_C = 40.0
- Cyclic deterioration parameter for accelerated reloading stiffness deterioration, lambda_A = 40.0
- Cyclic deterioration parameter for unloading stiffness deterioration, lambda_K = 15.0
- Rate of strength deterioration, c_S = 1.0
- Rate of post-capping strength deterioration, c_C = 1.0
- Rate of accelerated reloading stiffness deterioration, c_A = 1.0
- Rate of unloading stiffness deterioration, c_K = 1.0
- Pinching parameter: Ratio of force at break point to maximum force, kappa_F = 0.4
- Pinching parameter: Ratio of deformation at break point to residual plastic deformation, kappa_D = 0.2
- Rate of cyclic deterioration in the positive loading direction, D_pos = 1.0
- Rate of cyclic deterioration in the negative loading direction, D_neg = 1.0

- Assigned element mass, mass = 2.65e-5 

NOTE #1: For peak-oriented behavior, both kappa_F and kappa_D shall be taken as 1.0


++Example++

*Part, name="Stud Main"

*Node

      1,           0.,   -57.1500015,           0.
      
      2,           0.,   -57.1500015,           0.
      
*User element, nodes=2, type=VU1, properties=26, coordinates=3,

variables=38		

1, 2, 3	

2, 1, 2, 3	

*Element, type=VU1, elset=VUSPRING 		

1, 1, 2		

*UEL property, elset=VUSPRING

92000, 6.6, 11, 76000, 1.07, 0.2, 15.0, 10.2

5.0, 34000, 1.07, 0.2, 15.0, 40.0, 40.0, 40.0

15, 1.0, 1.0, 1.0, 1.0, 0.4, 0.2, 1.0

1.0, 2.65E-5

*End Part
**  


NOTE #2: A maximum of 8 UEL properties can be input per line

NOTE #3: The properties define the non-linear behavior of the spring in the loading direction. In the transverse
direction, the behavior is assumed to be elastic with a stiffness equal to Ke.


----------
AUTHOR																					   
----------

Code Developed by: Hammad El Jisr, École Polytechnique Fédérale de Lausanne (EPFL)


-----------
LICENSE																					   
-----------

This project is licensed under the MIT License - see the LICENSE.md file for details


--------------
REFERENCES																							   
--------------

[1] Ibarra, L. F., R. A. Medina, and H. Krawinkler. 2005. “Hysteretic models that incorporate strength and 
stiffness deterioration.” Earthquake Eng. Struct. Dyn. 34 (12): 1489–1511. https://doi.org/10.1002/eqe.495.

[2] El Jisr, H., Lignos, D. G., and Elkady, A. 2020. “Hysteretic Behavior of Moment-Resisting Frames 
Considering Slab Restraint and Framing Action.” Journal of Structural Engineering, 146(8).
https://doi.org/10.1061/(ASCE)ST.1943-541X.0002696


