// minimal DAE example (dense)
// next step is to add function to different file

// #include "residual.h" // add the residual header file
// #include "jacobian.h" // add the jacobian header file
// #include "events.h" // add the events header file

#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunmatrix/sunmatrix_sparse.h>       /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>           /* defs. of SUNRabs, SUNRexp, etc.      */

#define NEQ 2

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
namespace py = pybind11;

using residual_type = std::function<py::array_t<double>(double, py::array_t<double>, py::array_t<double>)>;
using jacobian_type = std::function<py::array_t<double>(double, py::array_t<double>, double)>;
using event_type = std::function<py::array_t<double>(double, py::array_t<double>)>;
using np_array = py::array_t<double>;

class PybammFunctions
{
public:
  PybammFunctions(const residual_type &res, const jacobian_type &jac, const event_type &event, const int n_s, int n_e)
      : py_res(res), py_jac(jac), py_event(event), number_of_states(n_s), number_of_events(n_e)
  {
  }

  int number_of_states;
  int number_of_events;

  py::array_t<double> operator()(double t, py::array_t<double> y, py::array_t<double> yp)
  {
    return py_res(t, y, yp);
  }

  py::array_t<double> res(double t, py::array_t<double> y, py::array_t<double> yp)
  {
    return py_res(t, y, yp);
  }

  py::array_t<double> jac(double t, py::array_t<double> y, double cj)
  {
    return py_jac(t, y, cj);
  }

  np_array events(double t, np_array y)
  {
    return py_event(t, y);
  }

private:
  residual_type py_res;
  jacobian_type py_jac;
  event_type py_event;
};

int residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
  PybammFunctions *python_functions_ptr = static_cast<PybammFunctions *>(user_data);
  PybammFunctions python_functions = *python_functions_ptr;

  realtype *yval, *ypval, *rval;
  yval = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);
  rval = N_VGetArrayPointer(rr);

  int n = python_functions.number_of_states;
  py::array_t<double> y_np = py::array_t<double>(n, yval);
  py::array_t<double> yp_np = py::array_t<double>(n, ypval);

  py::array_t<double> r_np;

  r_np = python_functions.res(tres, y_np, yp_np);

  double *r_np_ptr = (double *)r_np.request().ptr;

  // just copying data
  int i;
  for (i = 0; i < n; i++)
  {
    rval[i] = r_np_ptr[i];
  }
  return 0;
}

int jacobian(realtype tt, realtype cj,
             N_Vector yy, N_Vector yp, N_Vector resvec,
             SUNMatrix JJ, void *user_data,
             N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{

  realtype *yval, *ypval;
  yval = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);

  PybammFunctions *python_functions_ptr = static_cast<PybammFunctions *>(user_data);
  PybammFunctions python_functions = *python_functions_ptr;

  int n = python_functions.number_of_states;
  py::array_t<double> y_np = py::array_t<double>(n, yval);
  py::array_t<double> yp_np = py::array_t<double>(n, ypval);

  // create pointer to data
  realtype *jac_data_ptr = SM_DATA_D(JJ);

  py::array_t<double> jac_np_array;

  jac_np_array = python_functions.jac(tt, y_np, cj);

  double *jac_np_data_ptr = (double *)jac_np_array.request().ptr;

  // just copying data into Sunmatrix (figure out how to pass pointers later)
  int j_size = n * n;
  int i;
  for (i = 0; i < j_size; i++)
  {
    jac_data_ptr[i] = jac_np_data_ptr[i];
  }

  // SM_ELEMENT_D(JJ, 0, 0) = -cj;

  // sparse stuff
  // realtype *csr_index_vals_ptr = SM_INDEXVALS_S(JJ);
  // realtype *csr_index_ptrs_ptr = SM_INDEXPTRS_S(JJ);
  // I think we only need to update the data values??
  // double a = 1.0;
  // csr_data_ptr = python_functions.jac_(a, y_np, yp_np);

  return (0);
}

int events(realtype t, N_Vector yy, N_Vector yp, realtype *events_ptr,
           void *user_data)
{

  realtype *yval, *ypval;
  yval = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);

  py::array_t<double> y_np = py::array_t<double>(2, yval);
  py::array_t<double> yp_np = py::array_t<double>(2, ypval);

  PybammFunctions *python_functions_ptr = static_cast<PybammFunctions *>(user_data);
  PybammFunctions python_functions = *python_functions_ptr;

  py::array_t<double> events_np_array;

  events_np_array = python_functions.events(t, y_np);

  double *events_np_data_ptr = (double *)events_np_array.request().ptr;

  // just copying data into Sunmatrix (figure out how to pass pointers later)
  int number_of_events = 2;
  int i;
  for (i = 0; i < number_of_events; i++)
  {
    events_ptr[i] = events_np_data_ptr[i];
  }

  // y = yval[0];
  // a = yval[1];
  // gout[0] = y - RCONST(1.0);
  // gout[1] = a - RCONST(0.0); // dummy condition to show can pass vectors

  return (0);
}

/* main program */

np_array solve(np_array t_np, np_array y0_np, np_array yp0_np,
               residual_type res, jacobian_type jac, event_type event,
               int number_of_events, int use_jacobian)
{
  auto t = t_np.unchecked<1>();
  // auto y0 = y0_np.unchecked<1>();
  auto yp0 = yp0_np.unchecked<1>();

  int number_of_states;
  number_of_states = y0_np.request().size;

  // create class
  // SetUserData pass pointer to my_rhs
  // within residual.c cast void pointer to PyBaMMRHS class
  // call my_rhs(args)

  auto y00 = y0_np.request();
  double *y0 = (double *)y00.ptr;

  void *ida_mem;          // pointer to memory
  N_Vector yy, yp, avtol; // y, y', and absolute tolerance
  realtype rtol, *yval, *ypval, *atval;
  int iout, retval, retvalr;
  int rootsfound[number_of_events];
  SUNMatrix J;
  SUNLinearSolver LS;
  SUNNonlinearSolver NLS;

  // allocate vectors
  yy = N_VNew_Serial(number_of_states);
  yp = N_VNew_Serial(number_of_states);
  avtol = N_VNew_Serial(number_of_states);

  // set initial value
  yval = N_VGetArrayPointer(yy);
  ypval = N_VGetArrayPointer(yp);
  // yval[0] = RCONST(y0[0]);
  // yval[1] = RCONST(y0[1]);
  int i;
  for (i = 0; i < number_of_states; i++)
  {
    yval[i] = y0[i];
    ypval[i] = yp0[i];
  }

  // ypval[0] = RCONST(yp0(0));
  // ypval[1] = RCONST(yp0(1));

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

  for (i = 0; i < number_of_states; i++)
  {
    atval[i] = RCONST(1.0e-8); // nb: this can be set differently for each state
  }

  retval = IDASVtolerances(ida_mem, rtol, avtol);

  // set events
  retval = IDARootInit(ida_mem, number_of_events, events);

  // set pybamm functions by passing pointer to it
  PybammFunctions pybamm_functions(res, jac, event, number_of_states, number_of_events);
  void *user_data = &pybamm_functions;
  IDASetUserData(ida_mem, user_data);

  // set linear solver
  J = SUNDenseMatrix(number_of_states, number_of_states); // jacobian type (must be dense for dense solvers, p183 of ida_guide.pdf)
  LS = SUNLinSol_Dense(yy, J);
  retval = IDASetLinearSolver(ida_mem, LS, J);

  // sparse stuff  (must use sparse solvers e.g. KLU or SuperLUMT, p183 of ida_guide.pdf)
  // J = SUNSparseMatrix(2, 2, 2, CSR_MAT); // template jacobian
  // // set the indices values and pts of the jacobian (leaving the values unset)
  // realtype *csr_index_vals_ptr = SM_INDEXVALS_S(JJ);
  // realtype *csr_index_ptrs_ptr = SM_INDEXPTRS_S(JJ);

  if (use_jacobian == 1)
  {
    printf("\nSetting jacobian \n");
    retval = IDASetJacFn(ida_mem, jacobian);
  }

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
  SUNMatDestroy(J);
  N_VDestroy(avtol);
  N_VDestroy(yp);

  printf("t=%f, y=%f, a=%f \n", tret, yval[0], yval[1]);

  return py::array_t<double>(number_of_states, yval);
}

class Dense
{
public:
  Dense(np_array t_np, np_array y0_np, np_array yp0_np)
      : number_of_timesteps(t_np.request().size), number_of_states(y0_np.request().size)
  {

    // convert numpy array ICs into N_Vectors
    auto tt = t_np.request();
    t = (double *)tt.ptr;

    auto y0 = y0_np.unchecked<1>();
    auto yp0 = yp0_np.unchecked<1>();

    yy = N_VNew_Serial(number_of_states);
    yp = N_VNew_Serial(number_of_states);

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    int i;
    for (i = 0; i < number_of_states; i++)
    {
      yval[i] = y0[i];
      ypval[i] = yp0[i];
    }

    // allocate memory for solver
    ida_mem = IDACreate();

    // initialise solver
    retval = IDAInit(ida_mem, residual, t[0], yy, yp);
  }

  void link_python_functions(const residual_type &res, const jacobian_type jac, const event_type evt, const int n_events)
  {
    // store python functions in user data
    PybammFunctions pybamm_functions(res, jac, evt, number_of_states, n_events);
    void *user_data = &pybamm_functions;
    IDASetUserData(ida_mem, user_data);
  }

  void set_linear_solver()
  {
    A = SUNDenseMatrix(number_of_states, number_of_states);
    LS = SUNLinSol_Dense(yy, A);
    retval = IDASetLinearSolver(ida_mem, LS, A);
  }

  void set_events()
  {
    const int n_events = number_of_events;
    retval = IDARootInit(ida_mem, n_events, events);
  }

  void set_jacobian()
  {
    retval = IDASetJacFn(ida_mem, jacobian);
  }

  void set_tolerances(const double rel_tol, const np_array abs_tol_py)
  {
    auto abs_tol_unchecked = abs_tol_py.unchecked<1>();

    N_Vector abs_tol;
    realtype *abs_tol_ptr;
    abs_tol_ptr = N_VGetArrayPointer(abs_tol);

    int i;
    for (i = 0; i < number_of_states; i++)
    {
      abs_tol_ptr[i] = abs_tol_unchecked[i];
    }

    retval = IDASVtolerances(ida_mem, rel_tol, abs_tol);
  }

  np_array solve()
  {
    realtype tout, tret;
    realtype t_final = t[number_of_timesteps - 1];
    while (tret < t_final)
    {
      // solver options:
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
    N_VDestroy(yp);
    // check that abs_tol etc are not unfreed

    printf("t=%f, y=%f, a=%f \n", tret, yval[0], yval[1]);

    return py::array_t<double>(number_of_states, yval);
  }

  int retval;

private:
  // python functions
  residual_type residual_py;
  jacobian_type jacobian_py;
  event_type events_py;

  // input data
  int number_of_timesteps;
  int number_of_events;
  int number_of_states;
  bool use_jacobian;

  // ida
  void *ida_mem;
  double *t;
  N_Vector yy, yp;
  realtype *yval, *ypval;
  SUNMatrix A;
  SUNLinearSolver LS;
};

PYBIND11_MODULE(sundials, m)
{
  m.doc() = "sundials solvers"; // optional module docstring

  m.def("solve", &solve, "The solve function",
        py::arg("t"), py::arg("y0"), py::arg("yp0"), py::arg("res"), py::arg("jac"),
        py::arg("events"), py::arg("number_of_events"),
        py::arg("use_jacobian"),
        py::return_value_policy::take_ownership);

  py::class_<PybammFunctions>(m, "pybamm_functions")
      .def(py::init<const residual_type &, const jacobian_type &, const event_type &, const int &, const int &>())
      .def("__call__", &PybammFunctions::operator())
      .def("res", &PybammFunctions::res);

  py::class_<Dense>(m, "Dense")
      .def(py::init<np_array &, np_array &, np_array &>())
      .def("link_python_functions", &Dense::link_python_functions)
      .def("set_linear_solver", &Dense::set_linear_solver)
      .def("set_events", &Dense::set_events)
      .def("set_jacobian", &Dense::set_jacobian)
      .def("set_tolerances", &Dense::set_tolerances)
      .def("solve", &Dense::solve)
      .def_readwrite("return_value", &Dense::retval);
}
