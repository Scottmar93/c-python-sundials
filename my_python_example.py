import os
import ctypes

# will probably need to have a cmake file instead
# for window. also using bash...
os.system("echo Cleaning previous build")
os.system("make clean")
os.system("echo Compiling .c files...")
os.system("make")

fun = ctypes.CDLL("my_simple_example.so")

t = fun.main()

print("The final time is", t)
