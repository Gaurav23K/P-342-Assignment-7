"""
Gaurav Kanu    [1811067]
Assignment: 7
Question : 1(b)-->  dy/dx = 6 - 2y/x
"""
import you_can

# Opening the data files :
a = open("data_Q1b_(0.5).txt", "w+")
b = open("data_Q1b_(0.1)).txt", "w+")
c = open("data_Q1b_(0.05)).txt", "w+")
d = open("data_Q1b_(0.02).txt", "w+")


def fun(x, y):
    return 6 - (2 * y / x)  # Defining the mathematical function


x0 = 3
y0 = 1

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
3              1              
3.5            3.666666667    
4.0            5.619047619    
4.5            7.214285714    
5.0            8.611111111        
.                   .
.                   .
.                   .

2) For 0.1 -->

X               Y               
3              1              
3.1            1.533333333    
3.2            2.034408602    
3.3            2.507258065    
3.4            2.95530303         
.                   .
.                   .
.                   .

3) For 0.05 -->

X               Y               
3              1              
3.05           1.266666667    
3.1            1.525136612    
3.15           1.775938657    
3.2            2.019559652        
.                   .
.                   .
.                   .

4) For 0.02 -->

X               Y               
3              1              
3.02           1.106666667    
3.04           1.21200883     
3.06           1.316061345    
3.08           1.418857929        
.                   .
.                   .
.                   .

'''