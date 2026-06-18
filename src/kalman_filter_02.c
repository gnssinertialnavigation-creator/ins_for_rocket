#include "kalman_filter_02.h"

enum kalman_update_static_type
{
  gravity_type = 0,
  omega_type
};

static void kalman_pk_update( const int nmeas,
                              const double k_gain[N_KALMAN_02_STATES][nmeas],
                              const double H[nmeas][N_KALMAN_02_STATES],
                              double pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES] );
static void kalman_gravity_and_omega_update( int mode_type, double imu_noise_std, double imu_bias_std,
                                             const double vector_body[3], const double meas[3],
                                             const double ins_qb2e[4],
                                             double pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES],
                                             double xk[N_KALMAN_02_STATES] );

//---------------------------------------------------------
static void kalman_pk_update( const int nmeas,
                              const double k_gain[N_KALMAN_02_STATES][nmeas],
                              const double H[nmeas][N_KALMAN_02_STATES],
                              double pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES] )
{
  int idx, jdx;
  double pk_tmp[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  memcpy( &pk_tmp[0][0], &pk[0][0], sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES );
  double I_KH[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  memset( &I_KH[0][0] , 0.0, ( sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES ) );
  matrix_product( N_KALMAN_02_STATES, nmeas, N_KALMAN_02_STATES, k_gain, H, I_KH );

  for( idx=0; idx<N_KALMAN_02_STATES; idx++ )
  {
    for( jdx=0; jdx<N_KALMAN_02_STATES; jdx++ )
    {
      if( idx == jdx )
        I_KH[idx][jdx] = 1.0 - I_KH[idx][jdx];
      else
        I_KH[idx][jdx] = -I_KH[idx][jdx];
    }
  }
  matrix_product( N_KALMAN_02_STATES, N_KALMAN_02_STATES, N_KALMAN_02_STATES, I_KH, pk_tmp, pk );
}

//---------------------------------------------------------
static void kalman_gravity_and_omega_update( int mode_type, double imu_noise_std, double imu_bias_std,
                                             const double vector_body[3], const double meas[3],
                                             const double ins_qb2e[4],
                                             double pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES],
                                             double xk[N_KALMAN_02_STATES] )
{
  int idx, jdx;
  double pk_tmp[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  memcpy( &pk_tmp[0][0], &pk[0][0], sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES );
  double pht[N_KALMAN_02_STATES][3];
  double H[3][N_KALMAN_02_STATES];
  double H_t[N_KALMAN_02_STATES][3];
  memset( &H[0][0], 0.0, ( sizeof( double ) * 3 * N_KALMAN_02_STATES ) );
  // -L
  double C_b2e[3][3], C_e2b[3][3];
  quaternions_2_matrix( ins_qb2e, C_b2e );
  matrix_t( 3, 3, C_b2e, C_e2b );
  H[0][6] = -2.0 * ( ins_qb2e[2] * vector_body[2] - ins_qb2e[3] * vector_body[1] );
  H[0][7] = -2.0 * ( ins_qb2e[2] * vector_body[1] + ins_qb2e[3] * vector_body[2] );
  H[0][8] = -2.0 * ( ins_qb2e[1] * vector_body[1] + ins_qb2e[0] * vector_body[2] - 2.0 * ins_qb2e[2] * vector_body[0] );
  H[0][9] = -2.0 * ( ins_qb2e[1] * vector_body[2] - 2.0 * ins_qb2e[3] * vector_body[0] - ins_qb2e[0] * vector_body[1] );
  H[1][6] = -2.0 * ( ins_qb2e[3] * vector_body[0] - ins_qb2e[1] * vector_body[2] );
  H[1][7] = -2.0 * ( ins_qb2e[2] * vector_body[0] - ins_qb2e[0] * vector_body[2] - 2.0 * ins_qb2e[1] * vector_body[1] );
  H[1][8] = -2.0 * ( ins_qb2e[1] * vector_body[0] + ins_qb2e[3] * vector_body[2] );
  H[1][9] = -2.0 * ( ins_qb2e[0] * vector_body[0] + ins_qb2e[2] * vector_body[2] - 2.0 * ins_qb2e[3] * vector_body[1] );
  H[2][6] = -2.0 * ( ins_qb2e[1] * vector_body[1] - ins_qb2e[2] * vector_body[0] );
  H[2][7] = -2.0 * ( ins_qb2e[3] * vector_body[0] + ins_qb2e[0] * vector_body[1] - 2.0 * ins_qb2e[1] * vector_body[2] );
  H[2][8] = -2.0 * ( ins_qb2e[3] * vector_body[1] - 2.0 * ins_qb2e[2] * vector_body[2] - ins_qb2e[0] * vector_body[0] );
  H[2][9] = -2.0 * ( ins_qb2e[1] * vector_body[0] + ins_qb2e[2] * vector_body[1] );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
    {
      if( mode_type == gravity_type )
        H[idx][jdx+10] = -C_b2e[idx][jdx];
      else
        H[idx][jdx+13] = -C_b2e[idx][jdx];
    }
  }
  double inv_3by3[3][3];
  matrix_t( 3, N_KALMAN_02_STATES, H, H_t );
  matrix_product( N_KALMAN_02_STATES, N_KALMAN_02_STATES, 3, pk_tmp, H_t, pht );
  matrix_product( 3, N_KALMAN_02_STATES, 3, H, pht, inv_3by3 );
  double Q[3][3], CQ[3][3];
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_noise_std + imu_bias_std ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      inv_3by3[idx][jdx] += Q[idx][jdx];
  }

  if( invert3x3( inv_3by3 ) )
  {
    double k_gain[N_KALMAN_02_STATES][3];
    matrix_product( N_KALMAN_02_STATES, 3, 3, pht, inv_3by3, k_gain );
    double dz[3] = { 0 };
    matrix_vector_product( 3, N_KALMAN_02_STATES, H, xk, dz );
    for( idx=0; idx<3; idx++ )
      dz[idx] = meas[idx] - dz[idx];
    for( idx=0; idx<N_KALMAN_02_STATES; idx++ )
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
void pk_propagator_02( double Pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES],
                       double imu_force_noise_std, double imu_gyro_noise_std,
                       double imu_force_bias_std, double imu_gyro_bias_std,
                       const double gyro_body[3], const double f_body[3],
                       const double p_ecef[3], const double qu[4], const double dt )
{
  int idx, jdx;
  double F[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  double F_t[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  double FP[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  memset( &F[0][0], 0.0, ( sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES ) );
  for( idx=0; idx<N_KALMAN_02_STATES; idx++ )
    F[idx][idx] = 1.0;

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

  // L
  F[3][6] = 2.0 * dt * ( qu[2] * f_body[2] - qu[3] * f_body[1] );
  F[3][7] = 2.0 * dt * ( qu[2] * f_body[1] + qu[3] * f_body[2] );
  F[3][8] = 2.0 * dt * ( qu[1] * f_body[1] + qu[0] * f_body[2] - 2.0 * qu[2] * f_body[0] );
  F[3][9] = 2.0 * dt * ( qu[1] * f_body[2] - 2.0 * qu[3] * f_body[0] - qu[0] * f_body[1] );
  F[4][6] = 2.0 * dt * ( qu[3] * f_body[0] - qu[1] * f_body[2] );
  F[4][7] = 2.0 * dt * ( qu[2] * f_body[0] - qu[0] * f_body[2] - 2.0 * qu[1] * f_body[1] );
  F[4][8] = 2.0 * dt * ( qu[1] * f_body[0] + qu[3] * f_body[2] );
  F[4][9] = 2.0 * dt * ( qu[0] * f_body[0] + qu[2] * f_body[2] - 2.0 * qu[3] * f_body[1] );
  F[5][6] = 2.0 * dt * ( qu[1] * f_body[1] - qu[2] * f_body[0] );
  F[5][7] = 2.0 * dt * ( qu[3] * f_body[0] + qu[0] * f_body[1] - 2.0 * qu[1] * f_body[2] );
  F[5][8] = 2.0 * dt * ( qu[3] * f_body[1] - 2.0 * qu[2] * f_body[2] - qu[0] * f_body[0] );
  F[5][9] = 2.0 * dt * ( qu[1] * f_body[0] + qu[2] * f_body[1] );

  // C
  double C_b2e[3][3], C_e2b[3][3];
  quaternions_2_matrix( qu, C_b2e );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      F[idx+3][jdx+10] += C_b2e[idx][jdx] * dt;
  }
  matrix_t( 3, 3, C_b2e, C_e2b );

  // U
  F[6][7] += 0.5 * dt * ( -gyro_body[0] );
  F[6][8] += 0.5 * dt * ( -gyro_body[1] );
  F[6][9] += 0.5 * dt * ( -gyro_body[2] + OMEGA_EARTH );
  F[7][6] += 0.5 * dt * (  gyro_body[0] );
  F[7][8] += 0.5 * dt * (  gyro_body[2] + OMEGA_EARTH );
  F[7][9] += 0.5 * dt * ( -gyro_body[1] );
  F[8][6] += 0.5 * dt * (  gyro_body[1] );
  F[8][7] += 0.5 * dt * ( -gyro_body[2] - OMEGA_EARTH );
  F[8][9] += 0.5 * dt * (  gyro_body[0] );
  F[9][6] += 0.5 * dt * (  gyro_body[2] - OMEGA_EARTH );
  F[9][7] += 0.5 * dt * (  gyro_body[1] );
  F[9][8] += 0.5 * dt * ( -gyro_body[0] );

  // qX
  double qX[4][3];
  double qX_t[3][4];
  qX[0][0] = 0.5 * ( -qu[1] ); qX[0][1] = 0.5 * ( -qu[2] ); qX[0][2] = 0.5 * ( -qu[3] );
  qX[1][0] = 0.5 * (  qu[0] ); qX[1][1] = 0.5 * ( -qu[3] ); qX[1][2] = 0.5 * (  qu[2] );
  qX[2][0] = 0.5 * (  qu[3] ); qX[2][1] = 0.5 * (  qu[0] ); qX[2][2] = 0.5 * ( -qu[1] );
  qX[3][0] = 0.5 * ( -qu[2] ); qX[3][1] = 0.5 * (  qu[1] ); qX[3][2] = 0.5 * (  qu[0] );
  matrix_t( 4, 3, qX, qX_t );
  for( idx=0; idx<4; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      F[idx+6][jdx+13] += qX[idx][jdx] * dt;
  }

  matrix_t( N_KALMAN_02_STATES, N_KALMAN_02_STATES, F, F_t );
  matrix_product( N_KALMAN_02_STATES, N_KALMAN_02_STATES, N_KALMAN_02_STATES, F, Pk, FP );
  matrix_product( N_KALMAN_02_STATES, N_KALMAN_02_STATES, N_KALMAN_02_STATES, FP, F_t, Pk );

  double Q[3][3], CQ[3][3], qXQ[4][3], Q_qu[4][4];
  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_force_noise_std * dt ), 2 );
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
  matrix_product( 4, 3, 3, qX, Q, qXQ );
  matrix_product( 4, 3, 4, qXQ, qX_t, Q_qu );
  for( idx=0; idx<4; idx++ )
  {
    for( jdx=0; jdx<4; jdx++ )
      Pk[idx+6][jdx+6] += Q_qu[idx][jdx];
  }

  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_force_bias_std * dt ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      Pk[idx+10][jdx+10] += Q[idx][jdx];
  }

  memset( &Q[0][0], 0.0, ( sizeof( double ) * 3 * 3 ) );
  for( idx=0; idx<3; idx++ )
    Q[idx][idx] = pow( ( imu_gyro_bias_std * dt ), 2 );
  matrix_product( 3, 3, 3, C_b2e, Q, CQ );
  matrix_product( 3, 3, 3, CQ, C_e2b, Q );
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
      Pk[idx+13][jdx+13] += Q[idx][jdx];
  }

  if( DEBUGER_ON && KALMAN_FILTER_02_DEBUGER_ON )
  {
    matrix_show( N_KALMAN_02_STATES, N_KALMAN_02_STATES, F );
    matrix_show( N_KALMAN_02_STATES, N_KALMAN_02_STATES, Pk );
  }
}

//---------------------------------------------------------
int kalman_update_02( const double time_last, int measured_attitude,
                      const int n_meas, double ins_pos[3], double ins_vel[3],
                      double ins_qb2e[4], double ins_bf[3], double ins_bg[3],
                      double pk[N_KALMAN_02_STATES][N_KALMAN_02_STATES],
                      const double meas[16], const double f_body[3],
                      const double gy_body[3], double imu_force_noise_std,
                      double imu_gyro_noise_std, double imu_force_bias_std,
                      double imu_gyro_bias_std )
{
  int idx, jdx;
  double pk_tmp[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
  double xk[N_KALMAN_02_STATES] = { 0.0 };
  for( idx=0; idx<n_meas; idx++ )
  {
    memcpy( &pk_tmp[0][0], &pk[0][0], sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES );
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
    double kal_gain[N_KALMAN_02_STATES];
    for( jdx=0; jdx<N_KALMAN_02_STATES; jdx++ )
      kal_gain[jdx] = -1.0 * pk_tmp[jdx][idx] * inv_hpht;
    double dz = meas[idx] - ( -1.0 * xk[idx] );
    for( jdx=0; jdx<N_KALMAN_02_STATES; jdx++ )
      xk[jdx] += kal_gain[jdx] * dz;

    double I_KH[N_KALMAN_02_STATES][N_KALMAN_02_STATES];
    memset( &I_KH[0][0] , 0.0, ( sizeof( double ) * N_KALMAN_02_STATES * N_KALMAN_02_STATES ) );
    for( jdx=0; jdx<N_KALMAN_02_STATES; jdx++ )
      I_KH[jdx][idx] = ( -1.0 * kal_gain[jdx] );

    int kdx;
    for( jdx=0; jdx<N_KALMAN_02_STATES; jdx++ )
    {
      for( kdx=0; kdx<N_KALMAN_02_STATES; kdx++ )
    {
        if( jdx == kdx )
          I_KH[jdx][kdx] = 1.0 - I_KH[jdx][kdx];
        else
          I_KH[jdx][kdx] = -I_KH[jdx][kdx];
      }
    }
    matrix_product( N_KALMAN_02_STATES, N_KALMAN_02_STATES, N_KALMAN_02_STATES, I_KH, pk_tmp, pk );
  }

  if( ( time_last < 0.0 ) && ( n_meas > 6 ) )
  {
    // phi
    if( measured_attitude )
    {
      double pht[N_KALMAN_02_STATES][4];
      double H[4][N_KALMAN_02_STATES];
      double H_t[N_KALMAN_02_STATES][4];
      memset( &H[0][0], 0.0, ( sizeof( double ) * 4 * N_KALMAN_02_STATES ) );
      for( idx=0; idx<4; idx++ )
        H[idx][idx+6] = 1;
      matrix_t( 4, N_KALMAN_02_STATES, H, H_t );
      matrix_product( N_KALMAN_02_STATES, N_KALMAN_02_STATES, 4, pk, H_t, pht );
      double inv_4by4[4][4];
      matrix_product( 4, N_KALMAN_02_STATES, 4, H, pht, inv_4by4 );
      for( idx=0; idx<4; idx++ )
        inv_4by4[idx][idx] += pow( 0.05, 2 );

      if( lu_matrix_inv( 4, inv_4by4 ) )
      {
        double k_gain[N_KALMAN_02_STATES][4];
        matrix_product( N_KALMAN_02_STATES, 4, 4, pht, inv_4by4, k_gain );
        double dz[4] = { 0 };
        matrix_vector_product( 4, N_KALMAN_02_STATES, H, xk, dz );
        for( idx=0; idx<4; idx++ )
          dz[idx] = meas[idx+6] - dz[idx];
        for( jdx=0; jdx<N_KALMAN_02_STATES; jdx++ )
        {
          double sum = 0;
          for( idx=0; idx<4; idx++ )
            sum += k_gain[jdx][idx] * dz[idx];

          xk[jdx] += sum;
        }
        kalman_pk_update( 4, k_gain, H, pk );
      }
    }
    // f_b
    kalman_gravity_and_omega_update( gravity_type,
                                     imu_force_noise_std, imu_force_bias_std,
                                     &f_body[0], &meas[10],
                                     ins_qb2e, pk, xk );
    // gy_b
    kalman_gravity_and_omega_update( omega_type,
                                     imu_gyro_noise_std, imu_gyro_bias_std,
                                     &gy_body[0], &meas[13],
                                     ins_qb2e, pk, xk );
  }

  for( idx=0; idx<3; idx++ )
  {
    ins_pos[idx] -= xk[idx];
    ins_vel[idx] -= xk[idx+3];
    ins_bf[idx] += xk[idx+10];
    ins_bg[idx] += xk[idx+13];
  }

  for( idx=0; idx<4; idx++ )
    ins_qb2e[idx] -= xk[idx+6];

  vector_norm( 4, &ins_qb2e[0], &ins_qb2e[0] );

  return 1;
}
