"""
Gaurav Kanu    [1811067]
Assignment: 7
Question : 3-->  d2y/dx2 = dy/dx + 1
"""
import you_can

file1 = open("data_Q3h_(0.01).txt", "w+")
file2 = open("data_Q3l_(0.01).txt", "w+")


def fun2(x, y, z):
    return z+1


def fun1(x, y, z):
    return z


x0 = 0
y0 = 1
zh0 = 2
zl0 = 0

xn = 1
yn = 2 * (2.71828 - 1)


h = 0.01

you_can.shooting_method(x0, y0, zh0, zl0, xn, yn, h, fun1, fun2, file1, file2)

'''
***************************** OUTPUT ************************************


















'''