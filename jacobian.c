#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                          /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>           /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h> /* access to Newton SUNNonlinearSolver  */
#include <sundials/sundials_types.h>          /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>           /* defs. of SUNRabs, SUNRexp, etc.      */

int jacobian(realtype tt, realtype cj,
             N_Vector yy, N_Vector yp, N_Vector resvec,
             SUNMatrix JJ, void *user_data,
             N_Vector tempv1, N_Vector tempv2, N_Vector tempv3)
{
    realtype *yval;

    yval = N_VGetArrayPointer(yy);

    SM_ELEMENT_D(JJ, 0, 0) = -cj;
    SM_ELEMENT_D(JJ, 0, 1) = RCONST(1.0);
    SM_ELEMENT_D(JJ, 1, 0) = RCONST(0.0);
    SM_ELEMENT_D(JJ, 1, 1) = RCONST(-1.0);

    return (0);
}
