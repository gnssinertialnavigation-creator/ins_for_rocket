#include "navigation_process_function.h"

static int days_in_month( int month, int year );
static double gmst_rad( const DateTime *utc_p );
static void euler_method_prediction( const int n, const double dt, const double states_0[n],
                                     const double states_dot[n], double states[n] );

//---------------------------------------------------------
static int days_in_month( int month, int year )
{
  int days[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
  if( month == 2 )
  {
    int leap = ( year % 4 == 0 && year % 100 != 0 ) || ( year % 400 == 0 );
    return leap ? 29 : 28;
  }
  return days[month-1];
}

//---------------------------------------------------------
static double gmst_rad( const DateTime *utc_p )
{
  double DUT1 = 0.1734; // from IERS.
  double ut1_hours = utc_p->hours + DUT1 / 3600.0;

  double JD = 367.0 * utc_p->year -
    floor( 7.0 * ( utc_p->year + floor( ( utc_p->month + 9.0 ) / 12.0 ) ) / 4.0 ) +
    floor( 275.0 * utc_p->month / 9.0 ) + utc_p->day + 1721013.5 + ut1_hours / 24.0;

  double T0 = ( JD - 2451545.0 ) / 36525.0;

  double theta_deg = 280.46061837 +
    360.98564736629 * ( JD - 2451545.0 ) +
    0.000387933 * T0 * T0 -
    T0 * T0 * T0 / 38710000.0;

  theta_deg = fmod( theta_deg, 360.0 );
  if( theta_deg < 0.0 )
    theta_deg += 360.0;

  return theta_deg * M_PI / 180.0;
}

//---------------------------------------------------------
void add_seconds_to_utc( DateTime *utc_p, const double seconds )
{
  utc_p->hours += seconds / 3600.0;

  int day_carry = (int)floor(utc_p->hours / 24.0);
  utc_p->hours -= day_carry * 24.0;

  if( utc_p->hours < 0.0 )
  {
    utc_p->hours += 24.0;
    day_carry--;
  }

  utc_p->day += day_carry;

  while( utc_p->day > days_in_month( utc_p->month, utc_p->year ) )
  {
    utc_p->day -= days_in_month( utc_p->month, utc_p->year );
    utc_p->month++;
    if( utc_p->month > 12 )
    {
      utc_p->month = 1;
      utc_p->year++;
    }
  }

  while( utc_p->day < 1 )
  {
    utc_p->month--;
    if( utc_p->month < 1 )
    {
      utc_p->month = 12;
      utc_p->year--;
    }
    utc_p->day += days_in_month( utc_p->month, utc_p->year );
  }
}

//---------------------------------------------------------
void simple_gravity_model( const double pos_ecef[3], double g_ecef[3] )
{
  int idx;
  double r_2 = 0;
  for( idx=0; idx<3; idx++ )
    r_2 += ( pos_ecef[idx] * pos_ecef[idx] );

  double r_mag = sqrt(r_2);
  double g_mag = GM_EARTH / r_2;

  // |g| = GM/r^2;
  for( idx=0; idx<3; idx++ )
    g_ecef[idx] = ( -1.0 * ( pos_ecef[idx] * g_mag ) ) / r_mag;
}

//---------------------------------------------------------
void ecef_states_dot( const int n_st, const double states[n_st],
                      const double gyro_body[3], const double f_body[3],
                      double states_dot[n_st] )
{
  int idx;
  double qu[4], p_ecef[3], v_ecef[3];
  double negative_2_omega_r_dot[3], omega_2_r[3];
  memcpy( &p_ecef[0], &states[0], sizeof( double ) * 3 );
  memcpy( &v_ecef[0], &states[3], sizeof( double ) * 3 );
  memcpy( &qu[0], &states[6], sizeof( double ) * 4 );
  negative_2_omega_r_dot[0] = ( 2.0* OMEGA_EARTH * v_ecef[1] );
  negative_2_omega_r_dot[1] = -( 2.0* OMEGA_EARTH * v_ecef[0] );
  negative_2_omega_r_dot[2] = 0;
  omega_2_r[0] = ( OMEGA_EARTH * OMEGA_EARTH * p_ecef[0] );
  omega_2_r[1] = ( OMEGA_EARTH * OMEGA_EARTH * p_ecef[1] );
  omega_2_r[2] = 0;

  double f_ecef[3], C_b2e[3][3], C_e2b[3][3];
  quaternions_2_matrix( &qu[0], C_b2e );
  matrix_t( 3, 3, C_b2e, C_e2b );
  matrix_vector_product( 3, 3, C_b2e, f_body, f_ecef );
  double g_ecef[3] = { 0, 0, 0 };
  simple_gravity_model( p_ecef, g_ecef );

  for( idx=0; idx<3; idx++ )
  {
    states_dot[idx] = v_ecef[idx];
    states_dot[idx+3] = ( g_ecef[idx] + f_ecef[idx] );
    states_dot[idx+3] += ( negative_2_omega_r_dot[idx] + omega_2_r[idx] );
  }

  states_dot[6] = 0.5 * ( ( -qu[1] ) * gyro_body[0] +
                          ( -qu[2] ) * gyro_body[1] +
                          ( -qu[3] ) * gyro_body[2] ) + 0.5 * OMEGA_EARTH * (  qu[3] );

  states_dot[7] = 0.5 * ( (  qu[0] ) * gyro_body[0] +
                          ( -qu[3] ) * gyro_body[1] +
                          (  qu[2] ) * gyro_body[2] ) + 0.5 * OMEGA_EARTH * (  qu[2] );

  states_dot[8] = 0.5 * ( (  qu[3] ) * gyro_body[0] +
                          (  qu[0] ) * gyro_body[1] +
                          ( -qu[1] ) * gyro_body[2] ) + 0.5 * OMEGA_EARTH * ( -qu[1] );

  states_dot[9] = 0.5 * ( ( -qu[2] ) * gyro_body[0] +
                          (  qu[1] ) * gyro_body[1] +
                          (  qu[0] ) * gyro_body[2] ) + 0.5 * OMEGA_EARTH * ( -qu[0] );
}

//---------------------------------------------------------
static void euler_method_prediction( const int n, const double dt, const double states_0[n],
                                     const double states_dot[n], double states[n] )
{
  int idx;
  for( idx=0;idx<n;idx++ )
    states[idx] = states_0[idx] + states_dot[idx] * dt;
}

//---------------------------------------------------------
void euler_states_prediction( const int n_st, const double dt, const double st0[n_st],
                              const double st_dot[n_st], double st[n_st] )
{
  euler_method_prediction( n_st, dt, st0, st_dot, st );
  vector_norm( 4, &st[6], &st[6] );
}

//---------------------------------------------------------
void rk4_states_propagator( const double time_last, double pos_ecef[3],
                            double vel_ecef[3], double qu[4],
                            const double f0_body[3], const double gyro0_body[3],
                            const double f_body[3], const double gyro_body[3],
                            const double dt )
{
  if( time_last < 0.0 )
    return;

  int idx;
  double st[N_NAV_STATES], st0[N_NAV_STATES], f05_body[3], gyro05_body[3];
  memcpy( &st0[0], &pos_ecef[0], sizeof( double ) * 3 );
  memcpy( &st0[3], &vel_ecef[0], sizeof( double ) * 3 );
  memcpy( &st0[6], &qu[0], sizeof( double ) * 4 );
  for( idx=0; idx<3; idx++ )
  {
    f05_body[idx] = 0.5 * f0_body[idx] + 0.5 * f_body[idx];
    gyro05_body[idx] = 0.5 * gyro0_body[idx] + 0.5 * gyro_body[idx];
  }

  double k1[N_NAV_STATES], k2[N_NAV_STATES], k3[N_NAV_STATES], k4[N_NAV_STATES];
  ecef_states_dot( N_NAV_STATES, st0, gyro0_body, f0_body, k1 );
  euler_states_prediction( N_NAV_STATES, ( 0.5 * dt ), st0, k1, st );
  ecef_states_dot( N_NAV_STATES, st, gyro05_body, f05_body, k2 );
  euler_states_prediction( N_NAV_STATES, ( 0.5 * dt ), st0, k2, st );
  ecef_states_dot( N_NAV_STATES, st, gyro05_body, f05_body, k3 );
  euler_states_prediction( N_NAV_STATES, dt, st0, k3, st );
  ecef_states_dot( N_NAV_STATES, st, gyro_body, f_body, k4 );

  double st_dot[N_NAV_STATES];
  for( idx=0; idx<N_NAV_STATES; idx++ )
    st_dot[idx] = ( 1.0 / 6.0 ) * ( k1[idx] + 2.0 * k2[idx] + 2.0 * k3[idx] + k4[idx] );

  euler_states_prediction( N_NAV_STATES, dt, st0, st_dot, st );

  memcpy( &pos_ecef[0], &st[0], sizeof( double ) * 3 );
  memcpy( &vel_ecef[0], &st[3], sizeof( double ) * 3 );
  memcpy( &qu[0], &st[6], sizeof( double ) * 4 );
}

//---------------------------------------------------------
void lla2ecef( const double lat, const double lon, const double h,
               double pos_ecef[3] )
{
  // uint of lat and lon are rad.
  double A_2 = 6378137.0 * 6378137.0;
  double B_2 = 6356752.3142 * 6356752.3142;
  double sin_lat = sin( lat );
  double cos_lat = cos( lat );
  double sin_lon = sin( lon );
  double cos_lon = cos( lon );
  double N = A_2 / sqrt( A_2 * cos_lat * cos_lat + B_2 * sin_lat * sin_lat );
  pos_ecef[0] = ( N + h ) * cos_lat * cos_lon;
  pos_ecef[1] = ( N + h ) * cos_lat * sin_lon;
  pos_ecef[2] = ( ( B_2 / A_2 ) * N + h ) * sin_lat;
}

//---------------------------------------------------------
void ecef2lla( const double pos_ecef[3],
               double* lat_p, double* lon_p, double* h_p )
{
  // uint of lat and lon are rad.
  double A = 6378137.0;
  double A_2 = 6378137.0 * 6378137.0;
  double B = 6356752.3142;
  double B_2 = 6356752.3142 * 6356752.3142;
  double E = ( A_2 - B_2 ) / A_2;
  double F = ( A_2 - B_2 ) / B_2;
  double P = sqrt( pow(pos_ecef[0], 2 ) + pow( pos_ecef[1], 2 ) );
  double T = atan2( ( pos_ecef[2] * A ) , ( P * B ) );
  double sin_T = sin( T );
  double cos_T = cos( T );
  ( *lat_p ) = atan2( ( pos_ecef[2] + F * B * sin_T * sin_T * sin_T ) ,
                     ( P - E * A * cos_T * cos_T * cos_T ) );
  while( ( *lat_p ) > ( 0.5 * M_PI ) )
    ( *lat_p ) -= ( 2 * M_PI );

  while( ( *lat_p ) <= ( -( 0.5 * M_PI ) ) )
    ( *lat_p ) += ( 2 * M_PI );

  double sin_lat = sin( ( *lat_p ) );
  double cos_lat = cos( ( *lat_p ) );
  double N = A_2 / sqrt( A_2 * cos_lat * cos_lat + B_2 * sin_lat * sin_lat );
  ( *h_p ) = ( P / cos_lat ) - N;

  ( *lon_p ) = atan2( pos_ecef[1] , pos_ecef[0] );

  while( ( *lon_p ) > M_PI )
    ( *lon_p ) -= ( 2 * M_PI );

  while( ( *lon_p ) <= ( -M_PI ) )
    ( *lon_p ) += ( 2 * M_PI );
}

//---------------------------------------------------------
void ecef2enu_mat( const double lat, const double lon,
                   double C_ecef_2_enu[3][3] )
{
  C_ecef_2_enu[0][0] = -sin( lon );
  C_ecef_2_enu[0][1] = cos( lon );
  C_ecef_2_enu[0][2] = 0;
  C_ecef_2_enu[1][0] = -cos( lon ) * sin( lat );
  C_ecef_2_enu[1][1] = -sin( lon ) * sin( lat );
  C_ecef_2_enu[1][2] = cos( lat );
  C_ecef_2_enu[2][0] = cos( lon ) * cos( lat );
  C_ecef_2_enu[2][1] = sin( lon ) * cos( lat );
  C_ecef_2_enu[2][2] = sin( lat );
}

//---------------------------------------------------------
void ecef2enu( const double lat, const double lon,
               double C_ecef_2_enu[3][3],
               const double ecef[3], double enu[3] )
{
  ecef2enu_mat( lat, lon, C_ecef_2_enu );
  matrix_vector_product( 3, 3, C_ecef_2_enu, ecef, enu );
}

//---------------------------------------------------------
void ecef_to_eci(const DateTime *utc_time_p,
                 const double p_ecef[3], const double v_ecef[3],
                 double p_eci[3], double v_eci[3])
{
  double gmst = gmst_rad( utc_time_p );
  double c = cos( gmst );
  double s = sin( gmst );

  p_eci[0] = c * p_ecef[0] - s * p_ecef[1];
  p_eci[1] = s * p_ecef[0] + c * p_ecef[1];
  p_eci[2] = p_ecef[2];

  double Rvx = c * v_ecef[0] - s * v_ecef[1];
  double Rvy = s * v_ecef[0] + c * v_ecef[1];

  v_eci[0] = Rvx - OMEGA_EARTH * p_eci[1];
  v_eci[1] = Rvy + OMEGA_EARTH * p_eci[0];
  v_eci[2] = v_ecef[2];
}
