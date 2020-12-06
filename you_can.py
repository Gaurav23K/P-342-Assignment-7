import math
import random
import time


# Defining Euler forward method :

def euler_forward(x0, y0, h, fun, n, file):
    file.write("X{:<15}Y{:<15}\n".format('', ''))
    file.write("{:<15}{:<15}\n".format(x0, y0))
    while round(x0,6) < n:
        y = y0 + h * fun(x0, y0)
        x0 = x0 + h
        y0 = y

        file.write("{:<15.6}{:<15.10}\n".format(x0, y0))
    file.close()


# Defining Runge-Kutta method:

def runge_kutta(x0, y0, z0, h, fun1, fun2, n, file):
    file.write("X{:<15}Y{:<15}\n".format('', ''))
    file.write("{:<15}{:<15}\n".format(x0, y0))
    while round(x0,6) < n:
        k1y = h * fun1(x0, y0, z0)
        k1z = h * fun2(x0, y0, z0)
        k2y = h * fun1((x0 + h/2), (y0 + k1y/2), (z0 + k1z/2))
        k2z = h * fun2((x0 + h/2), (y0 + k1y/2), (z0 + k1z/2))
        k3y = h * fun1((x0 + h/2), (y0 + k2y/2), (z0 + k2z/2))
        k3z = h * fun2((x0 + h/2), (y0 + k2y/2), (z0 + k2z/2))
        k4y = h * fun1((x0 + h), (y0 + k3y), (z0 + k3z))
        k4z = h * fun2((x0 + h), (y0 + k3y), (z0 + k3z))
        y = y0 + (1/6 * (k1y + (2 * k2y) + (2 * k3y) + k4y))
        z0 = z0 + (1/6 * (k1z + (2 * k2z) + (2 * k3z) + k4z))
        x0 = x0 + h
        y0 = y
        file.write("{:<15.6}{:<15.10}\n".format(x0, y0))
    return y0


# Defining Shooting method :

def shooting_method(x0, y0, zh0, zl0, xn, yn, h, fun1, fun2, file1, file2):
    yh = runge_kutta(x0, y0, zh0, h, fun1, fun2, xn, file1)
    yl = runge_kutta(x0, y0, zl0, h, fun1, fun2, xn, file2)
    print("yh yl = 1=", yh, yl)
    print("yn = ", yn)
    if yh > yn > yl:
        while abs(yh - yn) > 0.001 or abs(yl - yn) > 0.001:
            if abs(yh - yn) > abs(yn - yl):
                zh0 = zl0 + (((zh0 - zl0)/(yh - yl)) * (yn - yl))
                print("zh0 =", zh0)

            elif abs(yh - yn) < abs(yn - yl):
                zl0 = zl0 + (((zh0 - zl0) / (yh - yl)) * (yn - yl))
                print("zl0 =",zl0)
            yh = runge_kutta(x0, y0, zh0, h, fun1, fun2, xn, file1)
            yl = runge_kutta(x0, y0, zl0, h, fun1, fun2, xn, file2)
            print("yh yl = 2=", yh, yl)
    elif yh < yn and yl < yn:
        return print("please change your 'zh' ")

    elif yh > yn and yl > yn:
        return print("please change your 'zl' ")

    elif yh < yn < yl:
        return print("please change your 'zl' and 'zh'")

#######################################  OTHER FUNCTIONS  ###################################


# Defining Monte carlo integration :

def monte_carlo(a, b, f, n):
    X = []
    for i in range(n):
        r = random.random()
        r1 = a+((b-a)*r)
        X.append(r1)
    sum, sum1 = 0.0, 0.0

    for i in range(n):
        sum += f(X[i])
        sum1 += (f(X[i]))**2

    p = ((b-a)*sum)/n
    var = sum1/n-(sum/n)**2
    sd = math.sqrt(var)
    return p, sd

# Defining simpson's method :
def simpson(a, b, f, n):
    h = (b-a)/(2*n)
    c = 0
    z = 0
    # Odd terms
    for i in range(1, n+1):
        x = h*(-1+2*i)+a
        c += 4*f(x)
        # Fourth order differentiation for error
        fd2x = d4x(f, x)
        if abs(fd2x) > z:
            z = abs(fd2x)
    # Even terms
    for i in range(1, n):
        x = 2*h*i+a
        c += 2*f(x)
        f4dx = d4x(f, x)
        if abs(f4dx) > z:
            z = abs(f4dx)
    # Integral
    s = (h/3)*(f(a)+f(b)+c)
    # max error
    error = (((b-a)**5)/(180*(n**4)))*z
    return s, error


# Fourth order differentiation
def d4x(f, x):
    h = 0.00001
    f4x = (d2x(f, x+h)+d2x(f, x-h)-2*d2x(f, x))/(h**2)
    return f4x


# Defining trapezoidal method :
def trapezoidal(a, b, f, n):
    h = (b-a)/n
    c = 0
    z = 0
    y = a
    k = y
    for i in range(1, n+1):
        xi = y + i*h
        c += (h/2)*(f(k) + f(xi))
        k = xi
        # Double differentiation for error
        fddx = d2x(f, xi)
        if abs(fddx) > z:
            z = abs(fddx)
    # Max error
    error = (((b - a)**3)/(12*(n**2)))*z
    return c, error


# Defining mid_point method :
def mid_point(a, b, f, n):
    # Dividing integral range into N equal parts
    h = (b-a)/n
    c = 0
    z = 0
    for i in range(0, n):
        x = ((a+i*h)+(a+(i+1)*h))/2
        c += h*f(x)
        # Double differentiation for error
        f2dx = d2x(f, x)
        if abs(f2dx) > z:
            z = abs(f2dx)
    # Max error
    error = (((b-a)**3)/(24*(n**2)))*z
    return c, error

# Defining double differentiation :
def d2x(f, x):
    h = 0.00001
    f2dx = (f(x+h)+f(x-h)-2*f(x))/(h**2)
    return f2dx


class color:

    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'


# Defining bracketing :
def bracketing(a, b, func):
    i = 0
    while i < 12:  # Limiting no.of iterations to 12
        if func(a) * func(b) < 0:
            print("The bracketing has been done properly it required " + str(i) + " iterations")
            break
        elif func(a) * func(b) == 0:
            if func(a) == 0:
                print("{} is the root of the nonlinear equation.".format(a))
            if func(b) == 0:
                print("{} is the root of the nonlinear equation.".format(b))
        elif func(a) * func(b) > 0:
            if abs(func(a)) < abs(func(b)):
                a = a - 1.5 * (b - a)
            else:
                b = b + 1.5 * (b - a)
        i += 1
    if i > 12:
        print("It took more than 12 iterations to bracket the root! Please choose new points.")
    else:
        return a, b


# Defining the bisection method
def bisection(a, b, func):
    file_1 = open('Q1_itr(bisection).txt', 'w')
    a, b = bracketing(a, b, func)  # Recalling the function
    i = 1
    print(color.BOLD + "Bisection method" + color.END)
    file_1.write(color.BOLD + "Bisection method" + color.END + "\n{:>3}       {:>10} ".format("Iterations",
                                                                                              "Absolute error") + "\n")
    while abs(a - b) > 0.000001 and i < 200:  # Limiting no.of iterations to 12
        c = (a + b) / 2
        if func(a) * func(c) < 0:
            file_1.write("{:>3}         {:>10.20f} ".format(i, abs(c - b)) + "\n")
            b = c
        elif func(a) * func(c) > 0:
            file_1.write("{:>3}         {:>10.20f} ".format(i, abs(c - a)) + "\n")
            a = c
        else:
            print("While iterating through the points one of the point has landed on the root of the equation.")
            print(a, b)
            break
        i += 1
    return "The root has been found to be at " + str(b) + " and it required " + str(i - 1) + " iterations" \
                                                                                             "\n_________________________________________________________________________________."

# Defining regular falsi method
def regular_falsi(a, b, func):
    file_2 = open('Q1_itr(regular_falsi).txt', 'w')
    bracketing(a, b, func)
    i = 1
    print(color.BOLD + "False Position method." + color.END)
    file_2.write(color.BOLD + "False Position method." + color.END + "\n{:>3}       {:>10} ".format("Iterations",
                                                                                                    "Absolute error") + "\n")
    c = 1
    c_n1 = a
    while abs(c - c_n1) > 0.000001 and i < 200:  # Limiting no.of iterations to 200
        c = b - ((b - a) * func(b)) / (func(b) - func(a))
        if func(a) * func(c) < 0:
            c_n1 = b
            file_2.write("{:>3}         {:>10.20f} ".format(i, abs(c - c_n1)) + "\n")
            b = c
        elif func(a) * func(c) > 0:
            c_n1 = a
            file_2.write("{:>3}         {:>10.20f} ".format(i, abs(c - c_n1)) + "\n")
            a = c
        else:
            if func(a) == 0:
                print(a, " is the root.")
            elif func(b) == 0:
                print(b, " is the root.")
        i += 1
    file_2.write("{:>3}         {:>10.20f} ".format(i, abs(c - c_n1)) + "\n")
    return "The root has been found to be at " + str(c) + " and it required " + str(i) + " iterations." \
            "\n_________________________________________________________________________________. "

# Defining newton raphson method
def newton_raphson(x_o, func):
    file_3 = open('Q1_itr(Newton_Raphson).txt', 'w')
    h = 0.01
    der = (func(x_o+h) - func(x_o-h))/(2*h)  # Derivative of the function
    x = x_o - (func(x_o) / der)
    copy = 1
    i = 1
    print(color.BOLD + "Newton-Raphson method." + color.END)
    file_3.write(color.BOLD + "Newton-Raphson method." + color.END + "\n{:>3}       {:>10} ".format("Iterations",
                                                                                                    "Absolute error") + "\n")
    while abs(x - copy) > 0.000001 and i < 200:  # Limiting no.of iterations to 200
        copy = x
        x_o = x
        der = (func(x_o + h) - func(x_o - h)) / (2 * h)
        x = x_o - (func(x_o) / der)
        file_3.write("{:>3}         {:>10.20f} ".format(i, abs(x - copy)) + "\n")
        i += 1
    file_3.write("{:>3}         {:>10.20f} ".format(i, abs(x - copy)) + "\n")
    return "The root has been found to be at " + str(x) + " and it required " + str(i) + " iterations."

# Defining the polynomial equation
def P(c, i, x):
    return ((c[0-i] * pow(x, 4))+(c[1-i] * pow(x, 3))+(c[2-i] * pow(x, 2))+(c[3-i] * pow(x, 1))+c[4-i])

# Defining first derivative
def P_1(P, c, i, x):
    h = 0.5
    s = (P(c, i, x+h) - P(c, i, x-h))/(2*h)
    return s

#Defining second derivative
def P_2(P, c, i, x):
    h = 0.5
    s = (P(c, i, x+h)+P(c, i, x-h)-2*P(c, i, x))/(2*h**2)
    return s


# Defining Laguerre method
def laguerre_Method(P, c, i, a0):
    e = len(c) - 1
    epsilon = math.pow(10, -4)
    if(P(c, i, a0) == 0):
        return a0
    else:
        for j in range(200):
            g = (P_1(P, c, i, a0)) / P(c, i, a0)
            h = (math.pow(g, 2)) - ((P_2(P, c, i, a0)) / P(c, i, a0))
            d1 = (g + math.sqrt(abs((e-i-1)*(( (e-i) * h ) - math.pow(g, 2)))))
            d2 = (g - math.sqrt(abs((e-i-1)*(( (e-i) * h) - math.pow(g, 2)))))
            if(abs(d1) > abs(d2)):
                a = (e-i) / d1
            else:
                a = (e-i) / d2
            a1 = a0 - a
            if(abs(a1 - a0) < epsilon):
                return a0
            a0 = a1


# Defining syntetic division method

def synthetic_Division(P, c, a0, r):
    for i in range(len(c) - 2):
        x1 = laguerre_Method(P, c, i, a0)
        for j in range(3 - i):
            c[j+1] = c[j+1] + (x1 * c[j])
        c[len(c) -1 - i] = 0
        r.append(x1)


#Reading a Matrix from a file
def read_Matrix(fil , A):
    file = open(fil , 'r')
    for line in file:
        ns = line.split()
        no = [float(n) for n in ns]
        A.append(no)
    file.close()

#To print the Matrix
def write_Matrix(x):
    for r in range(len(x)):
        print(x[r])

#Factorial Method
def factorial(num):
    fact = 1
    if num < 0:
        print("Sorry, factorial does not exist for negative numbers")
    elif num == 0:
        return 1
    else:
        for i in range(1, num + 1):
            fact = fact * i
        return fact

#Function for partial pivoting the Augmented Matrix / Only Matrix
def partial_pivot(a, b):
    n= len(a)
    counter = 0
    for r in range(0 , n):
        if abs(a[r][r]) == 0:
            for r1 in range(r+1 , n):
                if abs(a[r1][r]) > abs(a[r][r]):
                    counter = counter + 1
                    for x in range(0,n):
                        d1= a[r][x]
                        a[r][x] = a[r1][x]
                        a[r1][x] = d1
                    if(b!= 0):
                         d1 = b[r]
                         b[r] = b[r1]
                         b[r1] = d1
    return counter

#Multiplication of Matrices
def matrix_Multiplication(a, b):
    m = len(b[0])
    l=len(b)
    n = len(a)
    p2 = [[0 for y in range(m)] for x in range(n)]
    for x in range(n):
        for i in range(m):
            for y in range(l):
                p2[x][i] = p2[x][i] + (a[x][y] * b[y][i])
    return p2

#Gauss-Jordan Method
def gauss_Jordan(a , b):
    n = len(a)
    bn = len(b[0])
    for r in range(0 , n):
        partial_pivot(a , b)
        pivot = a[r][r]
        for c in range(r , n):
            a[r][c] = a[r][c]/pivot
        for c in range(0 , bn):
            b[r][c] = b[r][c]/pivot
        for r1 in range(0, n):
            if r1==r or a[r1][r] == 0:
                continue
            else:
                factor = a[r1][r]
                for c in range(r , n):
                    a[r1][c] = a[r1][c] - factor*a[r][c]
                for c in range(0 , bn):
                    b[r1][c] = b[r1][c] - factor*b[r][c]

# Forward- Backward Substitution
def forwardbackward_Substitution(a, b):
    m = len(b[0])
    n = len(a)
    # forward substitution
    y = [[0 for y in range(m)] for x in range(n)]
    for i in range(n):
        for j in range(m):
            s = 0
            for k in range(i):
                s = s + a[i][k] * y[k][j]
            y[i][j] = b[i][j] - s
    # backward substitution
    x = [[0 for y in range(m)] for x in range(n)]
    for i in range(n-1, -1, -1):
        for j in range(m):
            s = 0
            for k in range(i + 1, n):
                s = s + a[i][k] * x[k][j]
            x[i][j] = (y[i][j] - s) / a[i][i]

    return x

#L-U decomposition
def lu_Decomposition(a):
    n = len(a)
    for i in range(n):
        for j in range(n):
            s = 0
            if(i<=j):
                for k in range(i):
                    s = s + (a[i][k] * a[k][j])
                a[i][j] = a[i][j] - s
            else:
                for k in range(j):
                    s = s + (a[i][k] * a[k][j])
                a[i][j] = (a[i][j] - s) / a[j][j]

#Finding determinant of Upper Triangular Matrix
def uppertriangular_Determinant(a):
    n = len(a)
    p = 1
    for i in range(n-2):
        p = p * a[i][i]
    p = p * (a[n-2][n-2] * a[n-1][n-1])
    return p