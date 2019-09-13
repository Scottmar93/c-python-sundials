import sundials
import numpy as np

# times and initial conditions
t = np.linspace(0, 0.2, 100)
y0 = np.array([0.0, 1.0])
yp0 = np.array([1.0, 0.0])


# Standard pybamm functions
def rhs(t, y):
    r = np.zeros((y.size,))
    r[0] = 5 * y[1]
    r[1] = 1 - y[1]
    return r


def jac(t, y):
    J = np.zeros((2, 2))
    J[0][0] = 0.0
    J[0][1] = 1.0
    J[1][0] = 0.0
    J[1][1] = -1.0
    return J


num_of_events = 2


def events(t, y):
    g = np.zeros((num_of_events,))
    g[0] = y[0] - 1.0
    g[1] = y[1] - 0.0
    return g


# we need to process a little to get these into the required form
# for sundials

# the residual and jacobian residual are required by sundials
def res(t, y, yp):
    # must be of form r = f(t, y) - y'
    r = rhs(t, y)
    r[0] += -yp[0]
    return r


def jac_res(t, y, cj):
    # must be of form j_res = (dr/dy) - (cj) (dr/dy')
    # cj is just the input parameter
    # see p68 of the ida_guide.pdf for more details

    # jac = df/dy
    J = jac(t, y)

    # need to add on - cj * (dr/dy')
    # in our case this is just -cj for the (i, i) entry for rhs equations
    # and 0 otherwise
    # I think this will just be the - cj * mass_matrix?
    J[0][0] += -cj
    return J


use_jac = 1  # 1 to use, 0 to turn off

time = sundials.solve(t, y0, yp0, res, jac_res, events, num_of_events, use_jac)

print(time)

