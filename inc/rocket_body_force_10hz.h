#ifndef ROCKET_BODY_FORCE_10HZ_H_
#define ROCKET_BODY_FORCE_10HZ_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "basic_const.h"

#define IMU_FX_INPUT_BUFF_SIZE ( 500 * 10 )

//---------------------------------------------------------
void get_rocket_body_force( double t, double C_b2e[3][3],
                            double p_ecef[3], double f_body[3] );

#ifdef __cplusplus
}
#endif

#endif // #ifndef ROCKET_BODY_FORCE_10HZ_H_
