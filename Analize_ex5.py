# Ex5
# Aviva Malako 318995479
# Shir El Halimi 206085601
def inverse_matrix(matrix):
    """
    the methos  inverse mat using gauss elimination with
    :param matrix: list, mat
    :return: list, the inverse mat
    """
    temp_matrix = [[] for _ in matrix]
    for i, row in enumerate(matrix):
        assert len(row) == len(matrix)
        temp_matrix[i].extend(row + [0] * i + [1] + [0] * (len(matrix) - i - 1))
    gauss(temp_matrix)
    b = []
    for i in range(len(temp_matrix)):
        b.append(temp_matrix[i][len(temp_matrix[i]) // 2:])
    return b


def mul_matrix(mat_a, mat_b):
    """
    the method the calculate mul between the matrix
    :param mat_a: matrix
    :param mat_b:  matrix
    :return: the mul matrix
    """
    result_mulMatrix = []

    for i in range(0, len(mat_a)):
        tmp = []

        for j in range(0, len(mat_b[0])):
            sum = 0

            for k in range(0, len(mat_a[0])):
                sum += mat_a[i][k] * mat_b[k][j]
            tmp.append(sum)
        result_mulMatrix.append(tmp)

    return result_mulMatrix


def eliminate(x1, x2, col, target=0):
    """
    :param x1:  row
    :param x2:  row
    :param col:  col
    :param target: int
    """
    fac = (x2[col] - target) / x1[col]
    for i in range(len(x2)):
        x2[i] -= fac * x1[i]


def gauss(mat):
    """
    :param mat: matrix
    :return matrix: matrix
    """
    for i in range(len(mat)):
        if mat[i][i] == 0:
            for j in range(i + 1, len(mat)):
                if mat[i][j] != 0:
                    mat[i], mat[j] = mat[j], mat[i]
                    break
            else:
                raise ValueError("Matrix is not invertible")
        for j in range(i + 1, len(mat)):
            eliminate(mat[i], mat[j], i)
    for i in range(len(mat) - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            eliminate(mat[i], mat[j], i)
    for i in range(len(mat)):
        eliminate(mat[i], mat[i], i, target=1)
    return mat


def find_x_y(point_x, table):
    """
    Finding the closer numbers next to point x

    :param point_x:  the point we want to find
    :param table table that include x and y values
    :return: the closest x to point_x ansd their y according to the table
    """
    # find x1 and x2
    for i in table[0]:
        if point_x > table[0][i]:
            x1 = table[0][i]
        elif point_x < table[0][i]:
            x2 = table[0][i]
            break

    # find y1 and y2
    y1 = table[1][x1]
    y2 = table[1][x2]

    return x1, x2, y1, y2


# 1
def Lineary(point_x, table):
    """
    Finding point value with lineary method
    :param point_x: the point we want to find
    :param table:table that include x and y values
    :return:  point using lineary method
    """
    x1, x2, y1, y2 = find_x_y(point_x, table)
    return round(((y1 - y2) / (x1 - x2)) * point_x + (y2 * x1 - y1 * x2) / (x1 - x2), 4)


# 2
def Polynomial(point_x, table):
    """
    Finding point value with Polynomial method
    :param point_x: the point we want to find
    :param table:table that include x and y values
    :return: result of the polynomail
    """
    size = len(table[0])
    # declare mat with all values is 1
    matrix_n = [[1 for i in range(size)] for j in range(size)]
    # enter mat new values according to the table
    i = 0
    for item in table[0]:
        while i < size:
            k = 1
            for j in range(1, size):
                matrix_n[i][j] = item ** k
                k += 1
            break
        i += 1

    # transpose table[1] to be able to multiply it with mat
    temp_1 = []
    for item in table[1]:
        tmp = []
        tmp.append(item)
        temp_1.append(tmp)

    # saving the result of multiply matrices
    vector_a = mul_matrix(inverse_matrix(matrix_n), temp_1)

    result_p = 0
    for i in range(len(vector_a)):
        if i == 0:
            result_p += vector_a[i][0]
        else:
            result_p += vector_a[i][0] * (point_x ** i)

    return round(result_p, 4)


# 3
def Lagrange(point_x, table):
    """
    Finding point value with Lagrange method
    :param point_x: the point we want to find
    :param table: table that include x and y values
    :return: sum
    """

    list_n = []
    size = len(table[0])

    for i in range(size):
        temp_1 = []
        temp_1.append(table[0][i])
        temp_1.append(table[1][i])
        list_n.append(temp_1)

    # lagrange method
    sum = 0
    result_l = []
    tmp = 1
    for i in range(len(list_n)):
        for j in range(len(list_n)):
            if i != j:
                tmp *= (point_x - list_n[j][0]) / (list_n[i][0] - list_n[j][0])
        result_l.append(tmp)
        tmp = 1
    for i in range(size):
        sum += result_l[i] * list_n[i][1]
    return sum


# 4
def Neville(point_x, table):
    """
    Finding point value with Lagrange method
    :param point_x:the point we want to find
    :param table :table that include x and y values
    :return:the polynomial of degree n
    """

    n = len(table[0])
    p = n * [0]
    for k in range(n):
        for i in range(n - k):
            if k == 0:
                p[i] = table[1][i]
            else:
                p[i] = ((point_x - table[0][i + k]) * p[i] + \
                        (table[0][i] - point_x) * p[i + 1]) / \
                       (table[0][i] - table[0][i + k])
    return p[0]


def main():
    # define  points
    point_x1 = 2.5
    point_x2 = 3
    point_x3 = 1.5

    # define tables from the lesson
    table_Linear = [[0, 1, 2, 3, 4, 5, 6], [0, 0.8415, 0.9093, 0.1411, -0.7568, -0.9589, -0.2794]]
    table_Polynomial = [[1, 2, 3], [0.8415, 0.9093, 0.1411]]
    table_Lagrange = [[1, 2, 4], [1, 0, 1.5]]
    table_Neville = [[1, 1.3, 1.6, 1.9, 2.2], [0.7651, 0.6200, 0.4554, 0.2818, 0.1103]]

    while True:
        print("Which method would you like to use\nEnter 1  for Linear\nEnter 2 for Polynomial\n"
              "Enter 3 for Lagrange\nEnter 4 for Neville\nEnter another number that not between 1-4 to EXIT\n")
        user_input = input()

        if user_input == "1":
            print("The result =", Lineary(point_x1, table_Linear))

        elif user_input == "2":
            print("The result=", Polynomial(point_x1, table_Polynomial))

        elif user_input == "3":
            print("The result=", Lagrange(point_x2, table_Lagrange))

        elif user_input == "4":
            result = Neville(point_x3, table_Neville)
            print("The result=", round(result, 4))

        else:
            print("EXIT!")
            break


main()
