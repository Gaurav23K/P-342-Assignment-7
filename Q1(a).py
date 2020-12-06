"""
Gaurav Kanu    [1811067]
Assignment: 7
Question : 1(a)-->  dy/dx = yln(y)/x
"""

import math
import you_can

# Opening the data files :
a = open("data_Q1a_(0.5).txt", "w+")
b = open("data_Q1a_(0.1)).txt", "w+")
c = open("data_Q1a_(0.05)).txt", "w+")
d = open("data_Q1a_(0.02).txt", "w+")


def fun(x, y):
    return (y * math.log(y)) / x  # Defining the mathematical function


x0 = 2
y0 = 2.71828

h1 = 0.5
h2 = 0.1
h3 = 0.05
h4 = 0.02
n = 50

you_can.euler_forward(x0, y0, h1, fun, n, a)
you_can.euler_forward(x0, y0, h2, fun, n, b)
you_can.euler_forward(x0, y0, h3, fun, n, c)
you_can.euler_forward(x0, y0, h4, fun, n, d)

'''
*********************************** OUTPUT *******************************

1) For 0.5 -->

X               Y               
2              2.71828        
2.5            3.397849543    
3.0            4.229060546    
3.5            5.245430579    
4.0            6.487366697    
.                   .
.                   .
.                   .

2) For 0.1 -->

X               Y               
2              2.71828        
2.1            2.854193909    
2.2            2.996739075    
2.3            3.146238857    
2.4            3.303032044    
.                   .
.                   .
.                   .

3) For 0.05 -->

X               Y               
2              2.71828        
2.05           2.786236954    
2.1            2.855871942    
2.15           2.927226393    
2.2            3.000342745    
.                   .
.                   .
.                   .

4) For 0.02 -->

X               Y               
2              2.71828        
2.02           2.745462782    
2.04           2.772916041    
2.06           2.800642469    
2.08           2.828644784    
.                   .
.                   .
.                   .

'''