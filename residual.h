#ifndef RESIDUAL_H_DEF
#define RESIDUAL_H_DEF

#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>  /* access to serial N_Vector            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */

int residual(realtype tres, N_Vector yy, N_Vector yp, N_Vector resval, void *user_data);

#endif