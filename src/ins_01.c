#include "ins_01.h"

static double pos_ecef[3];
static double vel_ecef[3];
static double qu_b2e[4];
static double bias_force[3];
static double bias_gyro[3];
static double Pk[N_KALMAN_STATES][N_KALMAN_STATES];

static void set_ins_pos( const double p_ecef[3] );
static void set_ins_vel( const double v_ecef[3] );
static void set_ins_qu_b2e( const double qu[4] );
static void set_ins_bias_force( const double bf[3] );
static void set_ins_bias_gyro( const double bg[3] );
static void reset_ins_pk( void );
static void ins_propagator( const double gy_b[3], const double f_b[3],
                            const double time_last, const double dt );

//---------------------------------------------------------
void get_ins_01_pos( double p_ecef[3] )
{
  memcpy( &p_ecef[0], &pos_ecef[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void get_ins_01_vel( double v_ecef[3] )
{
  memcpy( &v_ecef[0], &vel_ecef[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void get_ins_01_qu_b2e( double q[4] )
{
  memcpy( &q[0], &qu_b2e[0], sizeof( double ) * 4 );
}

//---------------------------------------------------------
void get_ins_01_bias_force( double bf[3] )
{
  memcpy( &bf[0], &bias_force[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void get_ins_01_bias_gyro( double bg[3] )
{
  memcpy( &bg[0], &bias_gyro[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void reset_ins_01( const double p_ecef[3], const double v_ecef[3], const double qu[4] )
{
  set_ins_pos( p_ecef );
  set_ins_vel( v_ecef );
  set_ins_qu_b2e( qu );
  memset( &bias_force[0], 0.0, sizeof( double ) * 3 );
  memset( &bias_gyro[0], 0.0, sizeof( double ) * 3 );
  reset_ins_pk();
}

//---------------------------------------------------------
void ins_01( const double time_last, const double p_ecef[3], const double v_ecef[3],
             const double C_b2e_measured[3][3],
             const double f_body[3], const double gyro_body[3], const double dt )
{
  int idx;
  double gy_b[3], f_b[3], bf[3], bg[3], meas[15];
  get_ins_01_bias_force( &bf[0] );
  get_ins_01_bias_gyro( &bg[0] );
  for( idx=0; idx<3; idx++ )
  {
    f_b[idx] = f_body[idx] - bf[idx];
    gy_b[idx] = gyro_body[idx] - bg[idx];
  }
  ins_propagator( gy_b, f_b, time_last, dt );

  double pos[3], vel[3], qu[4];
  get_ins_01_pos( &pos[0] );
  get_ins_01_vel( &vel[0] );
  get_ins_01_qu_b2e( &qu[0] );
  for( idx=0; idx<3; idx++ )
  {
    meas[idx] = p_ecef[idx] - pos[idx];
    meas[idx+3] = v_ecef[idx] - vel[idx];
  }
  int n_meas = 6;
  double f_ecef[3] = {0}, gy_ecef[3] = {0};
  int measured_attitude = 0;
  if( time_last < 0.0 )
  {
    double g_ecef[3] = { 0, 0, 0 };
    double omega_ecef[3] = { 0, 0, OMEGA_EARTH };
    simple_gravity_model( p_ecef, g_ecef );

    double C_b2e[3][3];
    quaternions_2_matrix( &qu[0], C_b2e );
    matrix_vector_product( 3, 3, C_b2e, f_b, f_ecef );
    matrix_vector_product( 3, 3, C_b2e, gy_b, gy_ecef );

    if( measured_attitude )
    {
      n_meas += 3;
      double C_e2b[3][3], C_error[3][3];
      matrix_t( 3, 3, C_b2e, C_e2b );
      // ( I - phx )
      matrix_product( 3, 3, 3, C_b2e_measured, C_e2b, C_error );
      meas[6] = C_error[1][2];
      meas[7] = -C_error[0][2];
      meas[8] = C_error[0][1];
    }

    n_meas += 3;
    for( idx=0; idx<3; idx++ )
      meas[idx+9] = -g_ecef[idx] - f_ecef[idx];

    n_meas += 3;
    for( idx=0; idx<3; idx++ )
      meas[idx+12] = omega_ecef[idx] - gy_ecef[idx];
  }
  if( kalman_update( time_last, measured_attitude,
                     n_meas, pos, vel, qu, bf, bg, Pk,
                     meas, f_ecef, gy_ecef, imu_force_noise_std_01,
                     imu_gyro_noise_std_01, imu_force_bias_std_01,
                     imu_gyro_bias_std_01 ) )
  {
    set_ins_pos( &pos[0] );
    set_ins_vel( &vel[0] );
    set_ins_qu_b2e( &qu[0] );
    set_ins_bias_force( &bf[0] );
    set_ins_bias_gyro( &bg[0] );
  }
}

//---------------------------------------------------------
static void set_ins_pos( const double p_ecef[3] )
{
  memcpy( &pos_ecef[0], &p_ecef[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
static void set_ins_vel( const double v_ecef[3] )
{
  memcpy( &vel_ecef[0], &v_ecef[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
static void set_ins_qu_b2e( const double qu[4] )
{
  memcpy( &qu_b2e[0], &qu[0], sizeof( double ) * 4 );
}

//---------------------------------------------------------
static void set_ins_bias_force( const double bf[3] )
{
  memcpy( &bias_force[0], &bf[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
static void set_ins_bias_gyro( const double bg[3] )
{
  memcpy( &bias_gyro[0], &bg[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
static void reset_ins_pk( void )
{
  int idx;
  memset( Pk, 0.0, ( sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES )  );
  for( idx=0; idx<3; idx++ )
  {
    Pk[idx][idx] = 10000;
    Pk[idx+3][idx+3] = 4;
    Pk[idx+6][idx+6] = pow( ( 20 * DEG_2_RAD ), 2 );
    Pk[idx+9][idx+9] = pow( ( 10 * imu_force_bias_std_01 ), 2 );
    Pk[idx+12][idx+12] = pow( ( 10 * imu_gyro_bias_std_01 ), 2 );
  }
}

//---------------------------------------------------------
static void ins_propagator( const double gy_b[3], const double f_b[3],
                            const double time_last, const double dt )
{
  double st[N_NAV_STATES], st0[N_NAV_STATES], st_dot[N_NAV_STATES];
  get_ins_01_pos( &st0[0] );
  get_ins_01_vel( &st0[3] );
  get_ins_01_qu_b2e( &st0[6] );

  ecef_states_dot( N_NAV_STATES, st0, gy_b, f_b, st_dot );
  euler_states_prediction( N_NAV_STATES, dt, st0, st_dot, st );

  set_ins_pos( &st[0] );
  set_ins_vel( &st[3] );
  set_ins_qu_b2e( &st[6] );

  pk_propagator( Pk, imu_force_noise_std_01, imu_gyro_noise_std_01,
                 imu_force_bias_std_01, imu_gyro_bias_std_01, gy_b,
                 f_b, &st0[0], &st0[6], dt );
}
