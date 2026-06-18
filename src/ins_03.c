#include "ins_03.h"

static double pos_ecef[3];
static double vel_ecef[3];
static double qu_b2e[4];
static double bias_force[3];
static double bias_gyro[3];
static double Pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES];

static void set_ins_pos( const double p_ecef[3] );
static void set_ins_vel( const double v_ecef[3] );
static void set_ins_qu_b2e( const double qu[4] );
static void set_ins_bias_force( const double bf[3] );
static void set_ins_bias_gyro( const double bg[3] );
static void reset_ins_pk( void );
static void ins_propagator( const double gy_b[3], const double f_b[3],
                            const double time_last, const double dt );

//---------------------------------------------------------
void get_ins_03_pos( double p_ecef[3] )
{
  memcpy( &p_ecef[0], &pos_ecef[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void get_ins_03_vel( double v_ecef[3] )
{
  memcpy( &v_ecef[0], &vel_ecef[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void get_ins_03_qu_b2e( double q[4] )
{
  memcpy( &q[0], &qu_b2e[0], sizeof( double ) * 4 );
}

//---------------------------------------------------------
void get_ins_03_bias_force( double bf[3] )
{
  memcpy( &bf[0], &bias_force[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void get_ins_03_bias_gyro( double bg[3] )
{
  memcpy( &bg[0], &bias_gyro[0], sizeof( double ) * 3 );
}

//---------------------------------------------------------
void reset_ins_03( const double p_ecef[3], const double v_ecef[3], const double qu[4] )
{
  set_ins_pos( p_ecef );
  set_ins_vel( v_ecef );
  set_ins_qu_b2e( qu );
  memset( &bias_force[0], 0.0, sizeof( double ) * 3 );
  memset( &bias_gyro[0], 0.0, sizeof( double ) * 3 );
  reset_ins_pk();
}

//---------------------------------------------------------
void ins_03( const double time_last, const double p_ecef[3], const double v_ecef[3],
             const double qu_b2e_measured[4],
             const double f_body[3], const double gyro_body[3], const double dt )
{
  int idx;
  double gy_b[3], f_b[3], bf[3], bg[3], meas[16];
  get_ins_03_bias_force( &bf[0] );
  get_ins_03_bias_gyro( &bg[0] );
  for( idx=0; idx<3; idx++ )
  {
    f_b[idx] = f_body[idx] - bf[idx];
    gy_b[idx] = gyro_body[idx] - bg[idx];
  }
  ins_propagator( gy_b, f_b, time_last, dt );

  double pos[3], vel[3], qu[4];
  get_ins_03_pos( &pos[0] );
  get_ins_03_vel( &vel[0] );
  get_ins_03_qu_b2e( &qu[0] );
  for( idx=0; idx<3; idx++ )
  {
    meas[idx] = p_ecef[idx] - pos[idx];
    meas[idx+3] = v_ecef[idx] - vel[idx];
  }
  int n_meas = 6;
  double f_ecef[3], gy_ecef[3];
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
      n_meas += 4;
      for( idx=0; idx<4; idx++ )
        meas[idx+6] = -( qu_b2e_measured[idx] - qu[idx] );
    }

    n_meas += 3;
    for( idx=0; idx<3; idx++ )
      meas[idx+10] = -g_ecef[idx] - f_ecef[idx];

    n_meas += 3;
    for( idx=0; idx<3; idx++ )
      meas[idx+13] = omega_ecef[idx] - gy_ecef[idx];
  }
  if( kalman_update_02( time_last, measured_attitude,
                        n_meas, pos, vel, qu, bf, bg, Pk,
                        meas, f_b, gy_b, imu_force_noise_std_03,
                        imu_gyro_noise_std_03, imu_force_bias_std_03,
                        imu_gyro_bias_std_03 ) )
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
  memset( Pk, 0.0, ( sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES )  );
  for( idx=0; idx<3; idx++ )
  {
    Pk[idx][idx] = 10000;
    Pk[idx+3][idx+3] = 4;
    Pk[idx+10][idx+10] = pow( ( 10 * imu_force_bias_std_03 ), 2 );
    Pk[idx+13][idx+13] = pow( ( 10 * imu_gyro_bias_std_03 ), 2 );
  }
  for( idx=0; idx<4; idx++ )
    Pk[idx+6][idx+6] = pow( 0.1, 2 );
}

//---------------------------------------------------------
static void ins_propagator( const double gy_b[3], const double f_b[3],
                            const double time_last, const double dt )
{
  double st[N_NAV_STATES], st0[N_NAV_STATES], st_dot[N_NAV_STATES];
  get_ins_03_pos( &st0[0] );
  get_ins_03_vel( &st0[3] );
  get_ins_03_qu_b2e( &st0[6] );

  ecef_states_dot( N_NAV_STATES, st0, gy_b, f_b, st_dot );
  euler_states_prediction( N_NAV_STATES, dt, st0, st_dot, st );

  set_ins_pos( &st[0] );
  set_ins_vel( &st[3] );
  set_ins_qu_b2e( &st[6] );

  pk_propagator_02( Pk, imu_force_noise_std_03, imu_gyro_noise_std_03,
                    imu_force_bias_std_03, imu_gyro_bias_std_03, gy_b,
                    f_b, &st0[0], &st0[6], dt );
}
