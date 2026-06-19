#ifndef INS_02_H_
#define INS_02_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <string.h>
#include "basic_const.h"
#include "matrix_proc.h"
#include "kalman_filter.h"
#include "navigation_process_function.h"

#define imu_force_noise_std_02 ( 1e-2 * G_2_N )
#define imu_gyro_noise_std_02 ( 1e-2 * DEG_2_RAD )
#define imu_force_bias_std_02 ( 1e-3 * G_2_N )
#define imu_gyro_bias_std_02 ( 1e-3 * DEG_2_RAD )

void get_ins_02_pos( double p_ecef[3] );
void get_ins_02_vel( double v_ecef[3] );
void get_ins_02_qu_b2e( double q[4] );
void get_ins_02_bias_force( double bf[3] );
void get_ins_02_bias_gyro( double bg[3] );
void reset_ins_02( const double p_ecef[3], const double v_ecef[3], const double qu[4] );
void ins_02( const double time_last, const double p_ecef[3], const double v_ecef[3],
             const double C_b2e_measured[3][3],
             const double f_body[3], const double gyro_body[3], const double dt );

#ifdef __cplusplus
}
#endif

#endif // #ifndef INS_02_H_
