#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>  /* access to serial N_Vector            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */

int events(realtype t, N_Vector yy, N_Vector yp, realtype *gout,
           void *user_data)
{
    realtype *yval, y, a;

    yval = N_VGetArrayPointer(yy);
    y = yval[0];
    a = yval[1];
    gout[0] = y - RCONST(1.0);
    gout[1] = a - RCONST(0.0); // dummy condition to show can pass vectors

    return (0);
}