**SEMI-IMPLICIT METHOD FOR PRESSURE-LINKED EQUATION (SIMPLE) ALGORITHM**

This is a program to solve 2D flow field of steady-state, incompressible flow using SIMPLE algorithm and Staggered grid.

**ACKNOWLEDGEMENT**

Big thanks to 
- Yohanes Bimo Dwianto (https://www.linkedin.com/in/yohanes-bimo-dwianto-9b8b85119/)
- Mohamad Arif Andira (https://www.linkedin.com/in/andiraarif/)
- Muhammad Badruz (https://www.linkedin.com/in/muhammad-badruz-375856193/)

that share their time to discuss and help for the development of the code.


**CASE DESCRIPTION**

The domain is a rectangular with an inlet and outlet. You can specify the dimension of its length and width also the inlet and outlet dimension and position.

               _____________________             ___
              |                     |             |
              |                     |             |
              |                     |             |
           >>||                     |             |
        Li >>||                     |             | Ly
           >>||                     |             |
              |                     ||>>          |
              |                     ||>> Lo       |
              |_____________________||>>         _|_

              |---------------------|
                         Lx

**ASSUMPTIONS**
1. Steady state
2. Incompressible
3. Viscous
4. Laminar flow
5. No gravity

**RESULTS**

***Residual History***
![RMSE](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/RMSE.png) 

***Staggered Grid Arrangement***
![Grid Arrangement](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/Grid%20Arrangement.png) 

***Contour of U-momentum***
![u-contour](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/U%20momentum%20contour.png) 

***Contour of V-momentum***
![v-contour](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/V%20momentum%20contour.png) 

***Velocity Vector***
![velocity-vector](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/Velocity%20vector.png) 

***Horizontal u-momentum comparison***
![u-horizontal-comparison](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/horizontal%20u%20momentum%20comparison.png) 

***Horizontal v-momentum comparison***
![v-horizontal-comparison](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/horizontal%20v%20momentum%20comparison.png) 

***Vertical u-momentum comparison***
![u-vertical-comparison](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/vertical%20u%20momentum%20comparison.png) 

***Vertical v-momentum comparison***
![v-vertical-comparison](https://github.com/cahyaamalinadhi/SIMPLEAlgorithm/blob/master/datas/vertical%20v%20momentum%20comparison.png) 

