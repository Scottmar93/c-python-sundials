SHELL = sh

prefix       = /home/scott/Projects/ida_test/sundials/instdir
exec_prefix  = /home/scott/Projects/ida_test/sundials/instdir
includedir   = /home/scott/Projects/ida_test/sundials/instdir/include
libdir       = /home/scott/Projects/ida_test/sundials/instdir/lib

CPP      = /usr/bin/cc
CPPFLAGS = -O3 -DNDEBUG
CC       = /usr/bin/cc
CFLAGS   = -O3 -DNDEBUG
LDFLAGS  = 
LIBS     =  -lm /usr/lib/x86_64-linux-gnu/librt.so

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
LIBRARIES = -lsundials_idas -lsundials_nvecserial ${LIBS}

# -----------------------------------------------------------------------------------------

EXAMPLES = my_simple_example
EXAMPLES_DEPENDENCIES = residual

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------
all: python_shared_example 

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

${OBJECTS}: ${OBJECTS_DEPENDENCIES}

clean:
	rm -f ${OBJECTS_DEPENDENCIES} 
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}


