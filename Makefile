SHELL = sh

prefix       = /home/scott/Projects/ida_test/sundials/instdir
exec_prefix  = /home/scott/Projects/ida_test/sundials/instdir
includedir   = /home/scott/Projects/ida_test/sundials/instdir/include
libdir       = /home/scott/Projects/ida_test/sundials/instdir/lib

CPP      = /usr/bin/c++
CPPFLAGS = -O3 -DNDEBUG
CC       = /usr/bin/c++
CFLAGS   = -O3 -DNDEBUG
LDFLAGS  = 
LIBS     =  -lm /usr/lib/x86_64-linux-gnu/librt.so -lblas

LINKFLAGS = -Wl,-rpath,/home/scott/Projects/ida_test/sundials/instdir/lib

# -----------------------------------------------------------------------------------------

LIBRARIES_LAPACK = -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband  
LINKFLAGS_LAPACK = ${LINKFLAGS}::

# INCLUDES_KLU  = 
# LIBRARIES_KLU = -lsundials_sunlinsolklu 
# LINKFLAGS_KLU = ${LINKFLAGS}:

# INCLUDES_SLUMT  = 
# LIBRARIES_SLUMT = -lsundials_sunlinsolsuperlumt  
# LINKFLAGS_SLUMT = ${LINKFLAGS}::

TMP_INCS  = ${includedir} 
INCLUDES  = $(addprefix -I, ${TMP_INCS})
# INCLUDES  = $(${TMP_INCS})
LIBRARIES = -lsundials_idas -lsundials_nvecserial ${LIBS}

# -----------------------------------------------------------------------------------------

EXAMPLES = my_simple_example
EXAMPLES_DEPENDENCIES = residual jacobian events

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------
all: just_res

my_simple_example: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

python_shared_example: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -FPIC -shared -o $${i}.so $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -FPIC -shared -o $${i}.so $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

pybind: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` -I/usr/include/python3.6m -o $${i}`python3-config --extension-suffix`${INCLUDES}$ ${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS}  -L${libdir} ${LIBRARIES} ${LINKFLAGS} "; \
	  ${CC} -O3 -Wall -shared -std=c++11 -fPIC -I/usr/include/python3.6m `python3 -m pybind11 --includes` ${INCLUDES}$ ${i}.o ${OBJECTS_DEPENDENCIES} ${LDFLAGS}  -L${libdir} ${LIBRARIES} ${LINKFLAGS} -o $${i}`python3-config --extension-suffix`;
	done

just_res: 
	# g++ -O3 -Wall -shared -I/home/scott/Projects/ida_test/sundials/instdir/include -o residual.so residual.c -L/home/scott/Projects/ida_test/sundials/instdir/lib -lsundials_idas -lsundials_nvecserial -lm /usr/lib/x86_64-linux-gnu/librt.so -lblas -Wl,-rpath,/home/scott/Projects/ida_test/sundials/instdir/lib;
	g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` -I/home/scott/Projects/ida_test/sundials/instdir/include -o residual`python3-config --extension-suffix` residual.c -L/home/scott/Projects/ida_test/sundials/instdir/lib -lsundials_idas -lsundials_nvecserial -lm /usr/lib/x86_64-linux-gnu/librt.so -lblas -Wl,-rpath,/home/scott/Projects/ida_test/sundials/instdir/lib;

${OBJECTS}: ${OBJECTS_DEPENDENCIES}

clean:
	rm -f ${OBJECTS_DEPENDENCIES} 
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}.so



# c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` example.cpp -o example`python3-config --extension-suffix`