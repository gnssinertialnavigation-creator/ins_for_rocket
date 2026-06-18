#ifndef ROCKET_BODY_OMEGA_10HZ_H_
#define ROCKET_BODY_OMEGA_10HZ_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "basic_const.h"

#define IMU_GY_INPUT_BUFF_SIZE ( 500 * 10 )

//---------------------------------------------------------
void get_rocket_body_omega( double t, double C_b2e[3][3],
                            double gyro_body[3] );

#ifdef __cplusplus
}
#endif

#endif // #ifndef ROCKET_BODY_OMEGA_10HZ_H_
