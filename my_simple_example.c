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

/* main program */

double myFun(double t_end, double y0[], double yp0[])
{
  void *ida_mem;          // pointer to memory
  N_Vector yy, yp, avtol; // y, y', and absolute tolerance
  realtype rtol, *yval, *ypval, *atval;
  realtype t0, t, tout1, tout, tret;
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
  ypval[0] = RCONST(yp0[0]);
  ypval[1] = RCONST(yp0[1]);

  // set times
  t0 = RCONST(0.0);
  tout1 = RCONST(1.0);
  tout = RCONST(t_end);

  // allocate memory for solver
  ida_mem = IDACreate();

  // initialise solver
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

  t = RCONST(0.0);
  while (tret < tout)
  {
    // IDA_ONE_STEP_TSTOP
    // IDA_NORMAL
    retval = IDASolve(ida_mem, tout, &tret, yy, yp, IDA_ONE_STEP);

    if (retval == IDA_ROOT_RETURN)
    {
      break;
    }

    printf("t=%f, y=%f, a=%f \n", tret, yval[0], yval[1]);

    t = tret;
  }

  /* Free memory */
  IDAFree(&ida_mem);
  SUNLinSolFree(LS);
  SUNMatDestroy(A);
  N_VDestroy(avtol);
  N_VDestroy(yy);
  N_VDestroy(yp);

  return t;
}
