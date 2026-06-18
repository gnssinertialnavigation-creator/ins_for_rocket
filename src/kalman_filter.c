#include "kalman_filter.h"

enum kalman_update_static_type
{
  gravity_type = 0,
  omega_type
};

static void kalman_pk_update( const int nmeas,
                              const double k_gain[N_KALMAN_STATES][nmeas],
                              const double H[nmeas][N_KALMAN_STATES],
                              double pk[N_KALMAN_STATES][N_KALMAN_STATES] );
static void kalman_gravity_and_omega_update( int mode_type, double imu_noise_std, double imu_bias_std,
                                             const double vector_ecef[3], const double meas[3],
                                             const double C_b2e[3][3], const double C_e2b[3][3],
                                             double pk[N_KALMAN_STATES][N_KALMAN_STATES],
                                             double xk[N_KALMAN_STATES] );

//---------------------------------------------------------
static void kalman_pk_update( const int nmeas,
                              const double k_gain[N_KALMAN_STATES][nmeas],
                              const double H[nmeas][N_KALMAN_STATES],
                              double pk[N_KALMAN_STATES][N_KALMAN_STATES] )
{
  int idx, jdx;
  double pk_tmp[N_KALMAN_STATES][N_KALMAN_STATES];
  memcpy( &pk_tmp[0][0], &pk[0][0], sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES );
  double I_KH[N_KALMAN_STATES][N_KALMAN_STATES];
  memset( &I_KH[0][0] , 0.0, ( sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES ) );
  matrix_product( N_KALMAN_STATES, nmeas, N_KALMAN_STATES, k_gain, H, I_KH );

  for( idx=0; idx<N_KALMAN_STATES; idx++ )
  {
    for( jdx=0; jdx<N_KALMAN_STATES; jdx++ )
    {
      if( idx == jdx )
        I_KH[idx][jdx] = 1.0 - I_KH[idx][jdx];
      else
        I_KH[idx][jdx] = -I_KH[idx][jdx];
    }
  }
  matrix_product( N_KALMAN_STATES, N_KALMAN_STATES, N_KALMAN_STATES, I_KH, pk_tmp, pk );
}

//---------------------------------------------------------
static void kalman_gravity_and_omega_update( int mode_type, double imu_noise_std, double imu_bias_std,
                                             const double vector_ecef[3], const double meas[3],
                                             const double C_b2e[3][3], const double C_e2b[3][3],
                                             double pk[N_KALMAN_STATES][N_KALMAN_STATES],
                                             double xk[N_KALMAN_STATES] )
{
  int idx, jdx;
  double pk_tmp[N_KALMAN_STATES][N_KALMAN_STATES];
  memcpy( &pk_tmp[0][0], &pk[0][0], sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES );
  double pht[N_KALMAN_STATES][3];
  double H[3][N_KALMAN_STATES];
  double H_t[N_KALMAN_STATES][3];
  memset( &H[0][0], 0.0, ( sizeof( double ) * 3 * N_KALMAN_STATES ) );
  double tmp_3by3[3][3];
  skew_sym_3by3( vector_ecef, tmp_3by3 );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
    {
      H[idx][jdx+6] = tmp_3by3[idx][jdx];
      if( mode_type == gravity_type )
        H[idx][jdx+9] = -C_b2e[idx][jdx];
      else
        H[idx][jdx+12] = -C_b2e[idx][jdx];
    }
  }
  matrix_t( 3, N_KALMAN_STATES, H, H_t );
  matrix_product( N_KALMAN_STATES, N_KALMAN_STATES, 3, pk_tmp, H_t, pht );
  matrix_product( 3, N_KALMAN_STATES, 3, H, pht, tmp_3by3 );
  double Q[3][3], CQ[3][3];
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_noise_std + imu_bias_std ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      tmp_3by3[idx][jdx] += Q[idx][jdx];
  }

  if( invert3x3( tmp_3by3 ) )
  {
    double k_gain[N_KALMAN_STATES][3];
    matrix_product( N_KALMAN_STATES, 3, 3, pht, tmp_3by3, k_gain );
    double dz[3] = { 0 };
    matrix_vector_product( 3, N_KALMAN_STATES, H, xk, dz );
    for( idx=0; idx<3; idx++ )
      dz[idx] = meas[idx] - dz[idx];
    for( idx=0; idx<N_KALMAN_STATES; idx++ )
    {
      double sum = 0;
      for( jdx=0; jdx<3; jdx++ )
        sum += k_gain[idx][jdx] * dz[jdx];

      xk[idx] += sum;
    }
    kalman_pk_update( 3, k_gain, H, pk );
  }
}

//---------------------------------------------------------
void pk_propagator( double Pk[N_KALMAN_STATES][N_KALMAN_STATES],
                    double imu_force_noise_std, double imu_gyro_noise_std,
                    double imu_force_bias_std, double imu_gyro_bias_std,
                    const double gyro_body[3], const double f_body[3],
                    const double p_ecef[3], const double qu[4], const double dt )
{
  int idx, jdx;
  double F[N_KALMAN_STATES][N_KALMAN_STATES];
  double F_t[N_KALMAN_STATES][N_KALMAN_STATES];
  double FP[N_KALMAN_STATES][N_KALMAN_STATES];
  memset( &F[0][0], 0.0, ( sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES ) );
  for( idx=0; idx<N_KALMAN_STATES; idx++ )
    F[idx][idx] = 1.0;

  double C_b2e[3][3], C_e2b[3][3];
  quaternions_2_matrix( qu, C_b2e );
  matrix_t( 3, 3, C_b2e, C_e2b );

  // I
  F[0][3] += 1.0 * dt;
  F[1][4] += 1.0 * dt;
  F[2][5] += 1.0 * dt;

  // E
  double px_2 = p_ecef[0] * p_ecef[0];
  double py_2 = p_ecef[1] * p_ecef[1];
  double pz_2 = p_ecef[2] * p_ecef[2];
  double pxy = p_ecef[0] * p_ecef[1];
  double pyz = p_ecef[1] * p_ecef[2];
  double pxz = p_ecef[0] * p_ecef[2];
  double R = sqrt( p_ecef[0] * p_ecef[0] + p_ecef[1] * p_ecef[1] + p_ecef[2] * p_ecef[2] );
  double R_3 = pow( R, 3 ); double R_5 = pow( R, 5 );
  double inv_R_3 = 1.0 / R_3;
  F[3][0] += dt * GM_EARTH * ( 3.0 * px_2 / R_5 - inv_R_3 ) + OMEGA_EARTH * OMEGA_EARTH;
  F[3][1] += dt * GM_EARTH * 3.0 * pxy / R_5;
  F[3][2] += dt * GM_EARTH * 3.0 * pxz / R_5;
  F[4][0] += dt * GM_EARTH * 3.0 * pxy / R_5;
  F[4][1] += dt * GM_EARTH * ( 3.0 * py_2 / R_5 - inv_R_3 ) + OMEGA_EARTH * OMEGA_EARTH;
  F[4][2] += dt * GM_EARTH * 3.0 * pyz / R_5;
  F[5][0] += dt * GM_EARTH * 3.0 * pxz / R_5;
  F[5][1] += dt * GM_EARTH * 3.0 * pyz / R_5;
  F[5][2] += dt * GM_EARTH * ( 3.0 * pz_2 / R_5 - inv_R_3 );

  // 2omega
  F[3][4] += dt * ( 2.0 * OMEGA_EARTH );
  F[4][3] += dt * ( -2.0 * OMEGA_EARTH );

  // omega
  F[6][7] += dt * ( OMEGA_EARTH );
  F[7][6] += dt * ( -OMEGA_EARTH );

  // T
  double f_ecef[3], skew[3][3];
  matrix_vector_product( 3, 3, C_b2e, f_body, f_ecef );
  skew_sym_3by3( f_ecef, skew );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      F[idx+3][jdx+6] -= dt * skew[idx][jdx];
  }

  // C
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
    {
      F[idx+3][jdx+9] += dt * C_b2e[idx][jdx];
      F[idx+6][jdx+12] += dt * C_b2e[idx][jdx];
    }
  }

  matrix_t( N_KALMAN_STATES, N_KALMAN_STATES, F, F_t );
  matrix_product( N_KALMAN_STATES, N_KALMAN_STATES, N_KALMAN_STATES, F, Pk, FP );
  matrix_product( N_KALMAN_STATES, N_KALMAN_STATES, N_KALMAN_STATES, FP, F_t, Pk );
  double Q[3][3], CQ[3][3];
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_force_noise_std  * dt ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      Pk[idx+3][jdx+3] += Q[idx][jdx];
  }
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_gyro_noise_std * dt ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      Pk[idx+6][jdx+6] += Q[idx][jdx];
  }
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_force_bias_std * dt ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      Pk[idx+9][jdx+9] += Q[idx][jdx];
  }
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_gyro_bias_std * dt ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      Pk[idx+12][jdx+12] += Q[idx][jdx];
  }

  if( DEBUGER_ON && KALMAN_FILTER_DEBUGER_ON )
  {
    //matrix_show( N_KALMAN_STATES, N_KALMAN_STATES, F );
    matrix_show( N_KALMAN_STATES, N_KALMAN_STATES, Pk );
  }
}

//---------------------------------------------------------
int kalman_update( const double time_last, int measured_attitude,
                   const int n_meas, double ins_pos[3], double ins_vel[3],
                   double ins_qb2e[4], double ins_bf[3], double ins_bg[3],
                   double pk[N_KALMAN_STATES][N_KALMAN_STATES],
                   const double meas[15], const double f_ecef[3],
                   const double gy_ecef[3], double imu_force_noise_std,
                   double imu_gyro_noise_std, double imu_force_bias_std,
                   double imu_gyro_bias_std )
{
  int idx, jdx, kdx;
  double pk_tmp[N_KALMAN_STATES][N_KALMAN_STATES];
  double xk[N_KALMAN_STATES] = { 0.0 };
  double I_KH[N_KALMAN_STATES][N_KALMAN_STATES];
  for( idx=0; idx<6; idx++ )
  {
    memcpy( &pk_tmp[0][0], &pk[0][0], sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES );
    double hpht = pk_tmp[idx][idx];
    if( time_last < 0.0 )
    {
      if( idx < 3 )
        hpht += 1.0;
      else
        hpht += 0.01;
    }
    else
    {
      if( idx < 3 )
        hpht += 900.0;
      else
        hpht += 1.0;
    }

    double inv_hpht =  1.0 / hpht;
    double kal_gain[N_KALMAN_STATES];
    for( jdx=0; jdx<N_KALMAN_STATES; jdx++ )
      kal_gain[jdx] = -1.0 * pk_tmp[jdx][idx] * inv_hpht;
    double dz = meas[idx] - ( -1.0 * xk[idx] );
    for( jdx=0; jdx<N_KALMAN_STATES; jdx++ )
      xk[jdx] += kal_gain[jdx] * dz;

    memset( &I_KH[0][0] , 0.0, ( sizeof( double ) * N_KALMAN_STATES * N_KALMAN_STATES ) );
    for( jdx=0; jdx<N_KALMAN_STATES; jdx++ )
      I_KH[jdx][idx] = ( -1.0 * kal_gain[jdx] );

    for( jdx=0; jdx<N_KALMAN_STATES; jdx++ )
    {
      for( kdx=0; kdx<N_KALMAN_STATES; kdx++ )
      {
        if( jdx == kdx )
          I_KH[jdx][kdx] = 1.0 - I_KH[jdx][kdx];
        else
          I_KH[jdx][kdx] = -I_KH[jdx][kdx];
      }
    }
    matrix_product( N_KALMAN_STATES, N_KALMAN_STATES, N_KALMAN_STATES, I_KH, pk_tmp, pk );
  }

  double C_b2e[3][3], C_e2b[3][3];
  quaternions_2_matrix( ins_qb2e, C_b2e );
  matrix_t( 3, 3, C_b2e, C_e2b );
  if( ( time_last < 0.0 ) && ( n_meas > 6 ) )
  {
    // phi
    if( measured_attitude )
    {
      double pht[N_KALMAN_STATES][3];
      double H[3][N_KALMAN_STATES];
      double H_t[N_KALMAN_STATES][3];
      memset( &H[0][0], 0.0, ( sizeof( double ) * 3 * N_KALMAN_STATES ) );
      for( idx=0; idx<3; idx++ )
        H[idx][idx+6] = 1;
      matrix_t( 3, N_KALMAN_STATES, H, H_t );
      matrix_product( N_KALMAN_STATES, N_KALMAN_STATES, 3, pk, H_t, pht );
      double inv_3by3[3][3];
      matrix_product( 3, N_KALMAN_STATES, 3, H, pht, inv_3by3 );
      for( idx=0; idx<3; idx++ )
        inv_3by3[idx][idx] += pow( ( 0.5 * DEG_2_RAD ), 2 );

      if( invert3x3( inv_3by3 ) )
      {
        double k_gain[N_KALMAN_STATES][3];
        matrix_product( N_KALMAN_STATES, 3, 3, pht, inv_3by3, k_gain );
        double dz[3] = { 0 };
        matrix_vector_product( 3, N_KALMAN_STATES, H, xk, dz );
        for( idx=0; idx<3; idx++ )
          dz[idx] = meas[idx+6] - dz[idx];
        for( jdx=0; jdx<N_KALMAN_STATES; jdx++ )
        {
          double sum = 0;
          for( idx=0; idx<3; idx++ )
            sum += k_gain[jdx][idx] * dz[idx];

          xk[jdx] += sum;
        }
        kalman_pk_update( 3, k_gain, H, pk );
      }
    }
    // f_b
    kalman_gravity_and_omega_update( gravity_type,
                                     imu_force_noise_std, imu_force_bias_std,
                                     &f_ecef[0], &meas[9],
                                     C_b2e, C_e2b, pk, xk );
    // gy_b
    kalman_gravity_and_omega_update( omega_type,
                                     imu_gyro_noise_std, imu_gyro_bias_std,
                                     &gy_ecef[0], &meas[12],
                                     C_b2e, C_e2b, pk, xk );
  }

  for( idx=0; idx<3; idx++ )
  {
    ins_pos[idx] -= xk[idx];
    ins_vel[idx] -= xk[idx+3];
    ins_bf[idx] += xk[idx+9];
    ins_bg[idx] += xk[idx+12];
  }

  double _C[3][3], Cr_T[3][3];
  memcpy( &_C[0][0], &C_b2e[0][0], sizeof( double ) * 3 * 3 );
  double cr = cos( xk[6] ); double sr = sin( xk[6] );
  double cp = cos( xk[7] ); double sp = sin( xk[7] );
  double cy = cos( xk[8] ); double sy = sin( xk[8] );
  Cr_T[0][0] = cp*cy;         Cr_T[0][1] = cp*sy;         Cr_T[0][2] =   -sp;
  Cr_T[1][0] = sr*sp*cy-cr*sy;Cr_T[1][1] = cr*cy+sr*sp*sy;Cr_T[1][2] = sr*cp;
  Cr_T[2][0] = sr*sy+cr*sp*cy;Cr_T[2][1] = cr*sp*sy-sr*cy;Cr_T[2][2] = cr*cp;
  matrix_product( 3, 3, 3, Cr_T, _C, C_b2e );
  matrix_2_quaternions( C_b2e, ins_qb2e );

  return 1;
}
