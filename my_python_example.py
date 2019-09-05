import os
import ctypes

import numpy as np
import numpy.ctypeslib as npct

# will probably need to have a cmake file instead
# for window. also using bash...
os.system("echo Cleaning previous build")
os.system("make clean")
os.system("echo Compiling .c files...")
os.system("make")

# fun = ctypes.CDLL("my_simple_example.so")  # load using straight ctypes
fun = npct.load_library("my_simple_example.so", ".")

# set argument types
array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags="CONTIGUOUS")
fun.myFun.argtypes = [ctypes.c_double, array_1d_double, array_1d_double]

# set return types
# fun.myFun.restype = ctypes.c_double

fun.myFun.restype = ctypes.c_int

x = np.array([8.0, 6.0])

t = 10
y0 = np.array([0.0, 1.0])
yp0 = np.array([1.0, 0.0])
out = fun.myFun(t, y0, yp0)

print("The final is", out)
