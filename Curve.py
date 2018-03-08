
import numpy as np
import matplotlib.pyplot as plt


def factorial(n):
    """returns the factorial of number n"""
    if n == 0:
        return 1
    else:
        return n * factorial(n - 1)


def choose(n, r):
    """returns result of n choose r (nCr)"""
    return int(factorial(n) / (factorial(r) * factorial(n - r)))


def b(i, n, t):
    """returns the Bernstein weightings (b) of the control points in curve definition"""
    return choose(n, i) * t ** i * (1 - t) ** (n - i)


def bezier(n, t, a):
    """returns numerical answer of the Bezier curve definition given a t value and sets of control points"""
    result = 0

    for i in range(n + 1):
        result += a[i, :] * b(i, n, t)
    return result


def plot(a):
    """plot the curve given the control points"""
    n = len(a) - 1  # degree of curve
    t = np.linspace(0, 1, 101)
    bezier_plot = np.zeros([101, 2])

    for i in range(101):
        bezier_plot[i, :] = bezier(n, t[i], a)

    plt.plot(bezier_plot[:, 0], bezier_plot[:, 1])
    plt.plot(a[:, 0], a[:, 1], "ko")
    plt.plot(a[:, 0], a[:, 1], "k-")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


def bezier_given(a):
    """function that class gives"""
    # find order of curve from number of control points
    n = np.shape(a)[0]-1
    # initialise arrays
    b = np.zeros([101, 2])
    terms = np.zeros([n+1, 2])
    # create an array of values for t from 0 to 1 in 101 steps
    t = np.linspace(0, 1, 101)
    # loop through all t values
    for i in range(0, 101):
        # calculate terms inside sum in equation 13
        for j in range(0, n + 1):
            # YOUR CODE HERE
            terms[j, :] = a[j, :] * choose(n, j) * t[i] ** j * (1 - t[i]) ** (n - j)
        # sum terms to find Bezier curve
        b[i, :] = sum(terms, 0)
    # plot Bezier
    plt.plot(b[:, 0], b[:, 1])
    # plot control points
    plt.plot(a[:, 0], a[:, 1], 'ko')
    # plot control polygon
    plt.plot(a[:, 0], a[:, 1], 'k')
    return b


def rational_bezier(a, z):
    """rational bezier curve which allow the curve to go closer to a point according to the weight z"""
    # find order of curve from number of control points
    n = np.shape(a)[0] - 1
    # initialise arrays
    b_rat = np.zeros([101, 2])
    numerator_terms = np.zeros([n+1, 2])
    denominator_terms = np.zeros([n+1, 2])
    # create an array of values for t from 0 to 1 in 101 steps
    t = np.linspace(0, 1, 101)
    # loop through all t values
    for i in range(0, 101):
        # YOUR CODE HERE
        for j in range(0, n + 1):
            numerator_terms[j, :] = z[j] * a[j, :] * choose(n, j) * t[i] ** j * (1 - t[i]) ** (n - j)
            denominator_terms[j, :] = z[j] * choose(n, j) * t[i] ** j * (1 - t[i]) ** (n - j)
        b[i, :] = sum(numerator_terms, 0) / sum(denominator_terms, 0)
    # plot rational Bezier
    plt.plot(b_rat[:, 0], b_rat[:, 1])
    # plot control points
    plt.plot(a[:, 0], a[:, 1], 'ko')
    # plot control polygon
    plt.plot(a[:, 0], a[:, 1], 'k')
    return b_rat
