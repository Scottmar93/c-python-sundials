// minimal DAE example (dense)
// next step is to add function to different file

#include "residual.h" // add the residual header file
#include "jacobian.h" // add the jacobian header file
#include "events.h"   // add the events header file

#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>           /* defs. of SUNRabs, SUNRexp, etc.      */

#define NEQ 2

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
namespace py = pybind11;

/* main program */

py::array_t<double> solve(py::array_t<double> t_np, py::array_t<double> y0_np, py::array_t<double> yp0_np)
{
  auto t = t_np.unchecked<1>();
  // auto y0 = y0_np.unchecked<1>();
  auto yp0 = yp0_np.unchecked<1>();

  auto y00 = y0_np.request();
  double *y0 = (double *)y00.ptr;

  printf("t size %d", y00.size);
  int number_of_equations = y00.size;

  void *ida_mem;          // pointer to memory
  N_Vector yy, yp, avtol; // y, y', and absolute tolerance
  realtype rtol, *yval, *ypval, *atval;
  int iout, retval, retvalr;
  int rootsfound[2];
  SUNMatrix A;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;

  // allocate vectors
  yy = N_VNew_Serial(NEQ);
  yp = N_VNew_Serial(NEQ);
  avtol = N_VNew_Serial(NEQ);

  // set initial value
  yval = N_VGetArrayPointer(yy);
  yval[0] = RCONST(y0[0]);
  yval[1] = RCONST(y0[1]);

  ypval = N_VGetArrayPointer(yp);
  ypval[0] = RCONST(yp0(0));
  ypval[1] = RCONST(yp0(1));

  // // set times
  // t0 = RCONST(0.0);
  // tout1 = RCONST(1.0);
  // tout = RCONST(t_end);

  // allocate memory for solver
  ida_mem = IDACreate();

  // initialise solver
  realtype t0 = RCONST(0.0);
  retval = IDAInit(ida_mem, residual, t0, yy, yp);

  // set tolerances
  rtol = RCONST(1.0e-4);
  atval = N_VGetArrayPointer(avtol);
  atval[0] = RCONST(1.0e-8);
  atval[1] = RCONST(1.0e-6);

  retval = IDASVtolerances(ida_mem, rtol, avtol);

  // set events
  int num_of_events = 2;
  retval = IDARootInit(ida_mem, num_of_events, events);

  // set linear solver
  A = SUNDenseMatrix(NEQ, NEQ);
  LS = SUNLinSol_Dense(yy, A);
  retval = IDASetLinearSolver(ida_mem, LS, A);

  retval = IDASetJacFn(ida_mem, jacobian);

  realtype tout, tret;
  realtype t_final = t(99);

  while (tret < t_final)
  {
    // IDA_ONE_STEP_TSTOP
    // IDA_NORMAL
    retval = IDASolve(ida_mem, t_final, &tret, yy, yp, IDA_ONE_STEP);

    if (retval == IDA_ROOT_RETURN)
    {
      break;
    }

    printf("t=%f, y=%f, a=%f \n", tret, yval[0], yval[1]);
  }

  /* Free memory */
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(avtol);
  N_VDestroy(yp);

  printf("t=%f, y=%f, a=%f \n", tret, yval[0], yval[1]);

  return py::array_t<double>(number_of_equations, yval);
}

class PybammRHS
{
public:
  using function_type = std::function<double(double, double, double)>;

  PybammRHS(const function_type &f)
      : m_f(f)
  {
  }

  double operator()(double x, double f, double t)
  {
    return m_f(x, f, t);
  }

private:
  function_type m_f;
};

PYBIND11_MODULE(sundials, m)
{
  m.doc() = "sundials solvers"; // optional module docstring

  m.def("solve", &solve, "The solve function",
        py::arg("t"), py::arg("y0"), py::arg("yp0"),
        py::return_value_policy::take_ownership);

  py::class_<PybammRHS>(m, "RHS")
      .def(py::init<const PybammRHS::function_type &>())
      .def("__call__", &PybammRHS::operator());
}