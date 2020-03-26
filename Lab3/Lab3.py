import random
import itertools
import numpy as np
import scipy
from scipy.stats import f
from prettytable import PrettyTable


def f_critical(prob, f1, f2):
    return scipy.stats.f.ppf(prob, f1, f2)


def t_critical(prob, df):
    return scipy.stats.t.ppf(prob, df)


def c_critical(prob, f1, f2):
    return 1 / (1 + (f2 - 1) / scipy.stats.f.ppf(1 - (1 - prob) / f2, f1, (f2 - 1) * f1))


x_bounds = [[1, 1],
            [-30, 0],
            [-15, 35],
            [-30, -25]]
factors = len(x_bounds)
experiments = 4
samples = 2
confidence_prob = 0.9

combinations = list(itertools.product([-1, 1], repeat=4))
xn = [combinations[8],
      combinations[11],
      combinations[13],
      combinations[14]]
x = [[min(x_bounds[j]) if xn[i][j] < 0 else max(x_bounds[j])
      for j in range(len(x_bounds))]
     for i in range(len(xn))]

y_bounds = [int(200 + np.mean([min(x_bounds[i]) for i in range(1, factors)])),
            int(200 + np.mean([max(x_bounds[i]) for i in range(1, factors)])), ]


def create_matrix(y):
    table = PrettyTable()
    table_head = ["Екперемент №"]
    for i in range(factors):
        table_head.append(f"x{i}")

    for i in range(samples):
        table_head.append(f"y{i + 1}")

    table.field_names = table_head

    for i in range(experiments):
        table.add_row([i + 1, *xn[i], *y[i]])
    return table


def cochrans_test(samples):
    while True:
        y = [
            [random.randint(min(y_bounds), max(y_bounds)) for i in range(samples)]
            for j in range(experiments)
        ]
        matrix = create_matrix(y)

        s2_y = [np.var(y[i]) for i in range(experiments)]
        stat_c = max(s2_y) / sum(s2_y)
        crit_c = c_critical(confidence_prob, samples - 1, experiments)

        print(matrix)
        print(f"Розрахункова статистика С: {round(stat_c, 3)}")
        print(
            f"Критичний C для вірогідності довіри до {confidence_prob}: {round(crit_c, 3)}"
        )

        if stat_c < crit_c:
            print("Варіанти рівні.")
            break

        print("Варіанти не рівні. Збільшення розміру вибірки ...")
        samples += 1

    my = [np.mean(y[i]) for i in range(len(y))]

    xn_col = np.array(list(zip(*xn)))
    beta = [np.mean(my * xn_col[i]) for i in range(experiments)]

    yn = [sum(beta * np.array(xn[i])) for i in range(experiments)]

    print(f"Середні: {[round(my[i], 3) for i in range(experiments)]}")
    print(f"Обчислена функції: {[round(yn[i], 3) for i in range(experiments)]}")

    delta_x = [abs(x_bounds[i][0] - x_bounds[i][1]) / 2 for i in range(len(x_bounds))]
    x0 = [(x_bounds[i][0] + x_bounds[i][1]) / 2 for i in range(len(x_bounds))]
    b = [beta[0] - sum(beta[i] * x0[i] / delta_x[i] for i in range(1, factors))]
    b.extend([beta[i] / delta_x[i] for i in range(1, factors)])
    return s2_y, b, my, beta


def student_test(s2_y, beta, b):
    s2_b = np.mean(s2_y)
    s_beta = np.sqrt(s2_b / samples / experiments)
    stat_t = [abs(beta[i]) / s_beta for i in range(factors)]

    crit_t = t_critical(confidence_prob, (samples - 1) * experiments)

    print(f"Обчислена t статистика: {[round(stat_t[i], 3) for i in range(len(stat_t))]}")
    print(f"Критичний t для вірогідності довіри {confidence_prob}: {round(crit_t, 3)}")

    significant_coeffs = 0
    for i in range(len(stat_t)):
        if stat_t[i] < crit_t:
            b[i] = 0
            significant_coeffs += 1

    print(f"Коефіцієнти регресії: {[round(b[i], 3) for i in range(len(b))]}")

    y_calc = [sum((b * np.array(x))[i]) for i in range(experiments)]
    print(f"Розрахункові значення моделі: {[round(y_calc[i], 3) for i in range(len(y_calc))]}")
    return significant_coeffs, s2_b, y_calc


def fisher_test(s2_b, significant_coeffs, y_calc, my):
    s2_adeq = (samples / (experiments - significant_coeffs) * sum(
        [(y_calc[i] - my[i]) ** 2 for i in range(experiments)]))

    stat_f = s2_adeq / s2_b
    crit_f = f_critical(confidence_prob, (samples - 1) * experiments, experiments - significant_coeffs)

    print(f"Розрахункова статистика F: {round(stat_f, 3)}")
    print(f"Критичний F для вірогідності довіри до {confidence_prob}: {round(crit_f, 3)}")
    return stat_f, crit_f


s2_y, b, my, beta = cochrans_test(samples)
significant_coeffs, s2_b, y_calc = student_test(s2_y, beta, b)
stat_f, crit_f = fisher_test(s2_b, significant_coeffs, y_calc, my)

if stat_f > crit_f:
    print("Модель неадекватна.")
else:
    print("Модель адекватна.")
