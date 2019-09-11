import sundials
import numpy as np

t = np.linspace(0, 0.2, 100)
y0 = np.array([0.0, 1.0])
yp0 = np.array([1.0, 0.0])

time = sundials.solve(t, y0, yp0)

print(time)


def rhs(t, y, yp):

    out = np.zeros((3,))

    out[0] = y[0]
    out[1] = yp[1]
    out[2] = t

    return out


rhs_function = sundials.RHS(rhs)

t = 1
y = [10, 5]
y0 = [100, 50]

print(rhs_function(t, y, y0))

