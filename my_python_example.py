import sundials
import numpy as np

t = np.linspace(0, 0.2, 100)
y0 = np.array([0.0, 1.0])
yp0 = np.array([1.0, 0.0])

time = sundials.solve(t, y0, yp0)

print(time)


def my_fun(t, y, z):
    return t + y + z


rhs_function = sundials.RHS(my_fun)

print(rhs_function(1, 2, 3))

