#include <stdio.h>
#include <math.h>
#include <cblas.h>

#include <ida/ida.h>                   /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

int residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
    realtype *yval, *ypval, *rval;

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);

    double *a = yval + 1;
    printf("variable yval is at address: %p\n", (void *)&yval);
    printf("variable a is at address: %p\n", (void *)&a);
    printf("variable yval has value: %f\n", *yval);
    printf("variable a is has value: %f\n", *a);

    double A[] = {0.0, 1.0, 0.0, -1.0};
    double M[] = {-1.0, 0.0, 0.0, 0.0};
    double b[] = {0.0, 1.0};

    // b = A * y + b
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2, 2, 1.0, A, 2, yval, 1, 1.0, b, 1);

    // M * y' + b
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 2, 2, 1.0, M, 2, ypval, 1, 1.0, b, 1);

    rval = b;

    return 0;
}

int add(int i, int j)
{
    return i + j;
}

#include <pybind11/pybind11.h>
namespace py = pybind11;

PYBIND11_MODULE(residual, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function which adds two numbers",
          py::arg("i"), py::arg("j"));

    m.attr("the_answer") = 42;
    py::object world = py::cast("World");
    m.attr("what") = world;
} // namespace py=pybind11PYBIND11_MODULE(residual,m)