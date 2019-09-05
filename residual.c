#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>  /* access to serial N_Vector            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */

int residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector rr, void *user_data)
{
    realtype *yval, *ypval, *rval;

    yval = N_VGetArrayPointer(yy);
    ypval = N_VGetArrayPointer(yp);
    rval = N_VGetArrayPointer(rr);

    rval[0] = yval[1] - ypval[0];
    rval[1] = 1 - yval[1];

    return 0;
}