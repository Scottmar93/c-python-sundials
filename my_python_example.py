import sundials
import numpy as np

t = np.linspace(0, 0.2, 100)
y0 = np.array([0.0, 1.0])
yp0 = np.array([1.0, 0.0])


# Standard pybamm functions
def rhs(t, y):
    r = np.zeros((2,))
    r[0] = 10 * y[1]
    r[1] = 1 - y[1]
    return r


def jac(t, y):
    J = np.zeros((2, 2))
    J[0][0] = 0.0
    J[0][1] = 1.0
    J[1][0] = 0.0
    J[1][1] = -1.0
    return J


def events(t, y):
    g = np.zeros((2,))
    g[0] = y[0] - 1.0
    g[1] = y[1] - 0.0
    return g


# function residual and jac res required by sundials
def res(t, y, yp):
    r = rhs(t, y)
    r[0] += -yp[0]
    return r


def jac_res(t, y, cj):
    J = jac(t, y)
    J[0][0] += -cj
    return J


time = sundials.solve(t, y0, yp0, res, jac_res, events)

print(time)

