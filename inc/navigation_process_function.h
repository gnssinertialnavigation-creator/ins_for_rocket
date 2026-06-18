#ifndef NAVIGATION_PROCESS_FUNCTION_H_
#define NAVIGATION_PROCESS_FUNCTION_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "basic_const.h"
#include "matrix_proc.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  int    year;
  int    month;
  int    day;
  double hours;
} DateTime;

//---------------------------------------------------------
void add_seconds_to_utc( DateTime *utc_p, const double seconds );
//---------------------------------------------------------
void simple_gravity_model( const double pos_ecef[3], double g_ecef[3] );
//---------------------------------------------------------
void ecef_states_dot( const int n_st, const double states[n_st],
                      const double gyro_body[3], const double f_body[3],
                      double states_dot[n_st] );
//---------------------------------------------------------
void euler_states_prediction( const int n_st, const double dt, const double st0[n_st],
                              const double st_dot[n_st], double st[n_st] );
//---------------------------------------------------------
void rk4_states_propagator( const double time_last, double pos_ecef[3],
                            double vel_ecef[3], double qu[4],
                            const double f0_body[3], const double gyro0_body[3],
                            const double f_body[3], const double gyro_body[3],
                            const double dt );
//---------------------------------------------------------
void lla2ecef( const double lat, const double lon, const double h,
               double pos_ecef[3] );
//---------------------------------------------------------
void ecef2lla( const double pos_ecef[3],
               double* lat_p, double* lon_p, double* h_p );
//---------------------------------------------------------
void ecef2enu_mat( const double lat, const double lon,
                   double C_ecef_2_enu[3][3] );
//---------------------------------------------------------
void ecef2enu( const double lat, const double lon,
               double C_ecef_2_enu[3][3],
               const double ecef[3], double enu[3] );
//---------------------------------------------------------
void ecef_to_eci( const DateTime* utc_time_p,
                  const double p_ecef[3], const double v_ecef[3],
                  double p_eci[3], double v_eci[3] );

#ifdef __cplusplus
}
#endif

#endif // #ifndef NAVIGATION_PROCESS_FUNCTION_H_
