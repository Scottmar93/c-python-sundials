#ifndef EVENTS_H_DEF
#define EVENTS_H_DEF

#include <stdio.h>
#include <math.h>

#include <ida/ida.h>                 /* prototypes for IDA fcts., consts.    */
#include <nvector/nvector_serial.h>  /* access to serial N_Vector            */
#include <sundials/sundials_types.h> /* defs. of realtype, sunindextype      */

int events(realtype t, N_Vector yy, N_Vector yp, realtype *gout,
           void *user_data);

#endif