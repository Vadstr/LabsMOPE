import colorama
from colorama import Fore, Style

colorama.init()

m = 0
d = 0
N = 15

x1_min = -10
x1_max = 10
x2_min = -10
x2_max = 2
x3_min = -10
x3_max = 4

x01 = (x1_max - x1_min) / 2
x02 = (x2_max - x2_min) / 2
x03 = (x3_max - x3_min) / 2
delta_x1 = x1_max - x01
delta_x2 = x2_max - x02
delta_x3 = x3_max - x03
y_min = 200 + int((x1_min + x2_min + x3_min) / 3)
y_max = 200 + int((x1_max + x2_max + x3_max) / 3)

#  Ввід значень
correct_input = False
while not correct_input:
    try:
        m = int(input("Введіть кількість повторень: "))
        p = float(input("Введіть довірчу імовірність: "))
        correct_input = True
    except ValueError:
        pass

print(Style.RESET_ALL)

correct = False
while not correct:
    try:
        in_matrix = [
            [-1, -1, -1, +1, +1, +1, -1, +1, +1, +1],
            [-1, -1, +1, +1, -1, -1, +1, +1, +1, +1],
            [-1, +1, -1, -1, +1, -1, +1, +1, +1, +1],
            [-1, +1, +1, -1, -1, +1, -1, +1, +1, +1],
            [+1, -1, -1, -1, -1, +1, +1, +1, +1, +1],
            [+1, -1, +1, -1, +1, -1, -1, +1, +1, +1],
            [+1, +1, -1, +1, -1, -1, -1, +1, +1, +1],
            [+1, +1, +1, +1, +1, +1, +1, +1, +1, +1],
            [-1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0, 0],
            [+1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0, 0],
            [0, -1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0],
            [0, +1.215, 0, 0, 0, 0, 0, 0, 1.4623, 0],
            [0, 0, -1.215, 0, 0, 0, 0, 0, 0, 1.4623],
            [0, 0, +1.215, 0, 0, 0, 0, 0, 0, 1.4623],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        ]

        def out_red(text):
            print("\033[31m {}".format(text))

        class Critical_values:
            @staticmethod
            def get_cohren_value(size_of_selections, qty_of_selections, significance):
                from _pydecimal import Decimal
                from scipy.stats import f
                size_of_selections += 1
                partResult1 = significance / (size_of_selections - 1)
                params = [partResult1, qty_of_selections, (size_of_selections - 1 - 1) * qty_of_selections]
                fisher = f.isf(*params)
                result = fisher / (fisher + (size_of_selections - 1 - 1))
                return Decimal(result).quantize(Decimal('.0001')).__float__()

            @staticmethod
            def get_student_value(f3, significance):
                from _pydecimal import Decimal
                from scipy.stats import t
                return Decimal(abs(t.ppf(significance / 2, f3))).quantize(Decimal('.0001')).__float__()

            @staticmethod
            def get_fisher_value(f3, f4, significance):
                from _pydecimal import Decimal
                from scipy.stats import f
                return Decimal(abs(f.isf(significance, f4, f3))).quantize(Decimal('.0001')).__float__()


        def sqrt(element):
            from math import sqrt
            return sqrt(element)


        def fab(element):
            from math import fabs
            return fabs(element)


        def x(l1, l2, l3):
            x_1 = l1 * delta_x1 + x01
            x_2 = l2 * delta_x2 + x02
            x_3 = l3 * delta_x3 + x03
            return [x_1, x_2, x_3]


        def generate_matrix():
            from random import randrange
            matrix_with_y = [[randrange(y_min, y_max) for y in range(m)] for x in range(N)]
            return matrix_with_y


        def find_average(lst, orientation):
            average = []
            if orientation == 1:
                for rows in range(len(lst)):
                    average.append(sum(lst[rows]) / len(lst[rows]))
            else:
                for column in range(len(lst[0])):
                    number_lst = []
                    for rows in range(len(lst)):
                        number_lst.append(lst[rows][column])
                    average.append(sum(number_lst) / len(number_lst))
            return average


        def a(first, second):
            need_a = 0
            for j in range(N):
                need_a += x_matrix[j][first - 1] * x_matrix[j][second - 1] / N
            return need_a


        def find_known(number):
            need_a = 0
            for j in range(N):
                need_a += average_y[j] * x_matrix[j][number - 1] / 15
            return need_a


        def solve(lst_1, lst_2):
            from numpy.linalg import solve
            solver = solve(lst_1, lst_2)
            return solver


        def check_result(b_lst, k):
            y_i = b_lst[0] + b_lst[1] * matrix[k][0] + b_lst[2] * matrix[k][1] + b_lst[3] * matrix[k][2] + \
                  b_lst[4] * matrix[k][3] + b_lst[5] * matrix[k][4] + b_lst[6] * matrix[k][5] + b_lst[7] * matrix[k][6] + \
                  b_lst[8] * matrix[k][7] + b_lst[9] * matrix[k][8] + b_lst[10] * matrix[k][9]
            return y_i


        def student_test(b_lst, number_x=10):
            b_dispersion = sqrt(dispersion_b2)
            for column in range(number_x):
                t_practice = 0
                t_theoretical = Critical_values.get_student_value(f3, q)
                for row in range(N):
                    if column == 0:
                        t_practice += average_y[row] / N
                    else:
                        t_practice += average_y[row] * in_matrix[row][column - 1]
                if (t_practice / b_dispersion) < t_theoretical:
                    b_lst[column] = 0
            return b_lst


        def fisher_test():
            dispersion_ad = 0
            f4 = N - d
            for row in range(len(average_y)):
                dispersion_ad += (m * (average_y[row] - check_result(student_lst, row))) / (N - d)
            F_practice = dispersion_ad / dispersion_b2
            F_theoretical = Critical_values.get_fisher_value(f3, f4, q)
            return F_practice < F_theoretical


        x_matrix = [[] for x in range(N)]
        for i in range(len(x_matrix)):
            if i < 8:
                x1 = x1_min if in_matrix[i][0] == -1 else x1_max
                x2 = x2_min if in_matrix[i][1] == -1 else x2_max
                x3 = x3_min if in_matrix[i][2] == -1 else x3_max
            else:
                x_lst = x(in_matrix[i][0], in_matrix[i][1], in_matrix[i][2])
                x1, x2, x3 = x_lst
            x_matrix[i] = [x1, x2, x3, x1 * x2, x1 * x3, x2 * x3, x1 * x2 * x3, x1 ** 2, x2 ** 2, x3 ** 2]

        matrix_y = generate_matrix()
        average_x = find_average(x_matrix, 0)
        average_y = find_average(matrix_y, 1)
        matrix = [(x_matrix[i] + matrix_y[i]) for i in range(N)]
        mx_i = average_x
        my = sum(average_y) / 15

        values_arr = [
            [1, mx_i[0], mx_i[1], mx_i[2], mx_i[3], mx_i[4], mx_i[5], mx_i[6], mx_i[7], mx_i[8], mx_i[9]],
            [mx_i[0], a(1, 1), a(1, 2), a(1, 3), a(1, 4), a(1, 5), a(1, 6), a(1, 7), a(1, 8), a(1, 9), a(1, 10)],
            [mx_i[1], a(2, 1), a(2, 2), a(2, 3), a(2, 4), a(2, 5), a(2, 6), a(2, 7), a(2, 8), a(2, 9), a(2, 10)],
            [mx_i[2], a(3, 1), a(3, 2), a(3, 3), a(3, 4), a(3, 5), a(3, 6), a(3, 7), a(3, 8), a(3, 9), a(3, 10)],
            [mx_i[3], a(4, 1), a(4, 2), a(4, 3), a(4, 4), a(4, 5), a(4, 6), a(4, 7), a(4, 8), a(4, 9), a(4, 10)],
            [mx_i[4], a(5, 1), a(5, 2), a(5, 3), a(5, 4), a(5, 5), a(5, 6), a(5, 7), a(5, 8), a(5, 9), a(5, 10)],
            [mx_i[5], a(6, 1), a(6, 2), a(6, 3), a(6, 4), a(6, 5), a(6, 6), a(6, 7), a(6, 8), a(6, 9), a(6, 10)],
            [mx_i[6], a(7, 1), a(7, 2), a(7, 3), a(7, 4), a(7, 5), a(7, 6), a(7, 7), a(7, 8), a(7, 9), a(7, 10)],
            [mx_i[7], a(8, 1), a(8, 2), a(8, 3), a(8, 4), a(8, 5), a(8, 6), a(8, 7), a(8, 8), a(8, 9), a(8, 10)],
            [mx_i[8], a(9, 1), a(9, 2), a(9, 3), a(9, 4), a(9, 5), a(9, 6), a(9, 7), a(9, 8), a(9, 9), a(9, 10)],
            [mx_i[9], a(10, 1), a(10, 2), a(10, 3), a(10, 4), a(10, 5), a(10, 6), a(10, 7), a(10, 8), a(10, 9), a(10, 10)]
        ]
        known = [my, find_known(1), find_known(2), find_known(3), find_known(4), find_known(5), find_known(6), find_known(7),
                 find_known(8), find_known(9), find_known(10)]

        beta = solve(values_arr, known)
        print( "\tОтримане рівняння регресії")
        print(Style.RESET_ALL)
        print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
              "+ {:.3f} * Х1Х2X3 + {:.3f} * X11^2 + {:.3f} * X22^2 + {:.3f} * X33^2 = ŷ\n\t Перевірка"
              .format(beta[0], beta[1], beta[2], beta[3], beta[4], beta[5], beta[6], beta[7], beta[8], beta[9], beta[10]))
        for i in range(N):
            print("ŷ{} = {:.3f} ≈ {:.3f}".format((i + 1), check_result(beta, i), average_y[i]))

        y_dispersion = [0.0 for x in range(N)]
        homogeneity = False
        while not homogeneity:
            y_dispersion = [0.0 for x in range(N)]
            for i in range(N):
                dispersion_i = 0
                for j in range(m):
                    dispersion_i += (matrix_y[i][j] - average_y[i]) ** 2
                y_dispersion.append(dispersion_i / (m - 1))

            f1 = m - 1
            f2 = N
            f3 = f1 * f2
            q = 1 - p

            Gp = max(y_dispersion) / sum(y_dispersion)

            print("\n")
            print("\tКритерій Кохрена")
            print(Style.RESET_ALL)
            Gt = Critical_values.get_cohren_value(f2, f1, q)
            if Gt > Gp or m >= 25:
                print("\t\tДисперсія однорідна при рівні значимості {:.2f}!".format(q))
                homogeneity = True
            else:
                print("\t\tДисперсія не однорідна при рівні значимості {:.2f}!".format(q))
                m += 1
            if m == 25:
                exit()
        print(Style.RESET_ALL)
        dispersion_b2 = sum(y_dispersion) / (N * N * m)
        student_lst = list(student_test(beta))
        print("\tОтримане рівняння регресії з урахуванням критерія Стьюдента")
        print(Style.RESET_ALL)
        print("{:.3f} + {:.3f} * X1 + {:.3f} * X2 + {:.3f} * X3 + {:.3f} * Х1X2 + {:.3f} * Х1X3 + {:.3f} * Х2X3"
              "+ {:.3f} * Х1Х2X3 + {:.3f} * X11^2 + {:.3f} * X22^2 + {:.3f} * X33^2 = ŷ\n\tПеревірка"
              .format(student_lst[0], student_lst[1], student_lst[2], student_lst[3], student_lst[4], student_lst[5],
                      student_lst[6], student_lst[7], student_lst[8], student_lst[9], student_lst[10]))
        for i in range(N):
            print("ŷ{} = {:.3f} ≈ {:.3f}".format((i + 1), check_result(student_lst, i), average_y[i]))

        print("\n")
        print("\tКритерій Фішера")
        print(Style.RESET_ALL)
        d = 11 - student_lst.count(0)
        if fisher_test():
            correct = True
            print("\t\tРівняння регресії адекватне стосовно оригіналу")
        else:
            print("\t\tРівняння регресії неадекватне стосовно оригіналу")
    except ValueError:
        pass