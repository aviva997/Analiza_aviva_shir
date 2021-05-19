# Aviva Malako ,318995479
# Shir Halimi 206085607
import sympy as sp
from sympy.utilities.lambdify import lambdify

# epsilon
eps = 0.00000001


# A function for finding the values that cause the function to change a mark
def Find_renge(f, start_point, end_point):
    """
        Finds the suspected points and intersecting a function on the axis
        @:param: polynomial
        @:param: start point
        @:param end point
        @:return:  rannge list
    """
    # We run for all the point between the start point and the last point with jumps of 0.1
    range_list = []

    # While there are still range between the start point to the end point
    while end_point >= start_point:
        if f(start_point) * f(start_point + 0.1) < 0:
            range_list.append(round(start_point, 3))
            range_list.append(round(start_point + 0.1, 3))
        start_point += 0.1

    # Returns two lists that hold the values when they are placed to change the function
    return range_list


def Find_rengeft(f_d, start_point, end_point):
    """
        Finds the suspected points and intersecting a function derivative on the axis
        @:param: Derived from the polynomial
        @:param: start point
        @:param: end point
        @:return: range list  Derived
    """
    # We run for all the point between the start point and the last point with jumps of 0.1
    range_listd = []
    while end_point >= start_point:
        if f_d(start_point) * f_d(start_point + 0.1) < 0:
            range_listd.append(round(start_point, 3))
            range_listd.append(round(start_point + 0.1, 3))
        start_point += 0.1

    # Returns two lists that hold the values when they are placed to change the function
    return range_listd


def Bisection_Method(f, f_d, start_point, end_point, eps):
    """
        The function finds us the existing intersection points on the axis with Bisection method
        @:param: polynomial
        @:param:Derived from the polynomial
        @:param:  start point
        @:param: end point
        @:param: eps
        @:return: result_polinom
    """

    r_l_f = Find_renge(f, start_point, end_point)  # renge list f
    r_l_fd = Find_rengeft(f_d, start_point, end_point)  # renge list f_d
    counter = 0
    result_polinom = []

    if r_l_f and r_l_fd:
        for i in range(0, len(r_l_f), 2):

            # Whilw loop that run until the difference between the end point and the start point is greater than epsilon
            while (r_l_f[i + 1] - r_l_f[i]) > eps:

                x_m = (r_l_f[i] + r_l_f[i + 1]) / 2

                counter += 1
                print("{0} - {1}".format(counter, round(x_m, 6)))  # Printing the iterations
                if f(r_l_f[i]) * f(x_m) > 0:
                    r_l_f[i] = x_m

                else:
                    r_l_f[i + 1] = x_m

            if round(x_m, 6) not in result_polinom:
                print("x: ", round(x_m, 6))
                result_polinom.append(round(x_m, 6))  # Save the results

        for i in range(0, len(r_l_fd), 2):

            # Whilw loop that run until the difference between the end point and the start point is greater than epsilon
            while (r_l_fd[i + 1] - r_l_fd[i]) > eps:

                x_m = (r_l_fd[i] + r_l_fd[i + 1]) / 2

                counter += 1
                print("{0} - {1}".format(counter, round(x_m, 6)))  # Printing the iterations

                if f(r_l_fd[i]) * f(x_m) > 0:
                    r_l_fd[i] = x_m
                else:
                    r_l_fd[i + 1] = x_m

            if round(x_m, 6) not in result_polinom:
                print("x: ", round(x_m, 6))
                result_polinom.append(round(x_m, 6))

    n = len(result_polinom)
    i = 0
    while i < n and n != 0:  # list no empty
        if not -0.1 < f(result_polinom[i]) < 0.1:
            del result_polinom[i]
            n -= 1
            i = 0
        else:
            i += 1

    return result_polinom


def Secant_Method(f, start_point, end_point, eps):
    """
        The function finds us the existing intersection points on the axis with Secant method
         @:param: polynomial
        @:param:Derived from the polynomial
        @:param:  start point
        @:param: end point
        @:param: eps
        @:return: result_polinom
    """

    result_list = []
    counter = 1

    range_list = Find_renge(f, start_point, end_point)

    for i in range(0, len(range_list), 2):
        x = range_list[i]
        x_next = range_list[i + 1]

        while abs(x_next - x) > eps:

            print("{0} - {1}".format(counter, round(x, 6)))  # Printing the iterations
            if round(f(x_next), 4) == 0:
                print("x : ", round(x_next, 6))
                result_list.append(round(x, 6))

            last_r = x
            x = x_next
            x_next = (last_r * f(x) - x * f(last_r)) / (f(x) - f(last_r))
            counter += 1
    return result_list


def Newton_Raphson(f, f_d, start_point, end_point, eps):
    """
        The function finds us the existing intersection point  on the axis with Newton method
        notice - only one point then we run again from the main function
        @:param: polynomial
        @:param:Derived from the polynomial
        @:param:  start point
        @:param: end point
        @:param: eps
        @:return: result_polinom
    """

    x = (end_point + start_point) / 2
    if f(x) == 0:
        print(" trying  closer range")
        return None

    else:
        x_n = x - (f(x) / f_d(x))

    counter = 1

    if round(f(start_point), 5) == 0:
        print("x: ", round(start_point, 4))
        return round(x, 5)

    elif round(f(end_point), 5) == 0:
        print("x: ", round(end_point, 4))
        return end_point

    if f(start_point) > 0 and f(end_point) < 0 or f(start_point) < 0 and f(end_point) > 0:
        # A while loop that checks when the range between the two numbers is less than eps
        while x_n - x < eps:
            # Do it if the condition is not met

            print("{0} - {1}".format(counter, round(x, 6)))  # Printing the iterations

            if round(f(x), 5) == 0:
                print("x: ", round(x, 6))
                return round(x, 6)

            counter += 1
            x = x_n
            x_n = x - (f(x) / f_d(x))
    else:
        print("The function does not converge")


# main function
def Main():
    x = sp.symbols('x')
    start_point = 0
    end_Point = 0

    while True:

        i = int(input(
            "Enter 1 - to use the Bisection method \nEnter 2 - To use Newton Raphson \nEnter 3 - to use the Secant method\n"))

        if i == 1:

            print("\n**Bisection Method**")
            f1 = 4 * x ** 3 + 3 * x ** 2 - 6 * x
            f1_d = f1.diff(x)

            f1 = lambdify(x, f1)
            f1_d = lambdify(x, f1_d)

            start_Point = -2.0
            end_Point = 2.0
            result_list = Bisection_Method(f1, f1_d, start_Point, end_Point, eps)
            if result_list:

                print("The results of the function are: ", result_list)

            else:
                print("No result found!")

        elif i == 2:

            print("\n**Newton Raphson**")

            f2 = x ** 3 - x - 1
            # A new variable that holds the derivative of the function
            f_d2 = f2.diff(x)

            # order to insert x we will do
            f2 = lambdify(x, f2)
            f_d2 = lambdify(x, f_d2)

            result_list = []

            start_Point = 0.0
            end_Point = 3.0

            range_list_f = Find_renge(f2, start_Point, end_Point)
            range_list_fT = Find_rengeft(f_d2, start_Point, end_Point)

            for i in range(0, len(range_list_f), 2):
                answer = Newton_Raphson(f2, f_d2, range_list_f[i], range_list_f[i + 1], eps)
                if answer not in result_list and answer is not None:
                    result_list.append(answer)

            for i in range(0, len(range_list_fT), 2):
                answer = Newton_Raphson(f2, f_d2, start_Point, end_Point, eps)
                if answer not in result_list and answer is not None:
                    result_list.append(answer)

            if result_list:
                print("The results of the function are: ", result_list)

            else:
                print("No result found!")

        elif i == 3:

            print("\n**Secant Method**")
            y = sp.sin
            f3 = x ** 2 - y(x)

            # order to insert x we will do
            f3 = lambdify(x, f3)

            start_Point = -1.0
            end_Point = 1.0

            result_list = Secant_Method(f3, start_Point, end_Point, eps)
            if result_list:
                print("The results of the function are: ", result_list)

            else:
                print("No result found!")

        else:
            print("Exit")

            # stop the loop
            break


# Call to main function
Main()
