#ifndef KALMAN_FILTER_02_H_
#define KALMAN_FILTER_02_H_

#define N_KALMAN_02_STATES (16)

#define KALMAN_FILTER_02_DEBUGER_ON (0)

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <string.h>
#include "basic_const.h"
#include "matrix_proc.h"
#include "navigation_process_function.h"

void pk_propagator_02( double Pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES],
                       double imu_force_noise_std, double imu_gyro_noise_std,
                       double imu_force_bias_std, double imu_gyro_bias_std,
                       const double gyro_body[3], const double f_body[3],
                       const double p_ecef[3], const double qu[4], const double dt );
int kalman_update_02( const double time_last, int measured_attitude,
                      const int n_meas, double ins_pos[3], double ins_vel[3],
                      double ins_qb2e[4], double ins_bf[3], double ins_bg[3],
                      double pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES],
                      const double meas[16], const double f_body[3],
                      const double gy_body[3], double imu_force_noise_std,
                      double imu_gyro_noise_std, double imu_force_bias_std,
                      double imu_gyro_bias_std );

#ifdef __cplusplus
}
#endif

#endif // #ifndef KALMAN_FILTER_02_H_
