#include "matrix_proc.h"

//---------------------------------------------------------
void matrix_show( const int m, const int n, const double A[m][n] )
{
  int idx, jdx;
  printf( "mat,\n" );
  for( idx=0; idx<m; idx++ )
  {
    for( jdx=0; jdx<n; jdx++ )
    {
      if( A[idx][jdx] >= 0 )
        printf( " %.10f,", A[idx][jdx] );
      else
        printf( "%.10f,", A[idx][jdx] );
    }

    printf( "\n" );
  }
}

//---------------------------------------------------------
double vector_abs( const int n, const double in[n] )
{
  int idx;
  double abs = 0.0;
  for( idx=0; idx<n; idx++ )
    abs += in[idx] * in[idx];

  abs = sqrt( abs );
  return abs;
}

//---------------------------------------------------------
void vector_norm( const int n, const double in[n], double out[n] )
{
  int idx;
  memcpy( out, in, sizeof( double ) * n );
  double abs = vector_abs( n, out );
  for( idx=0; idx<n; idx++ )
    out[idx] /= abs;
}

//---------------------------------------------------------
void vector_sum( const int n, const double in_1[n], const double in_2[n], double out[n] )
{
  int idx;
  for( idx=0; idx<n; idx++ )
    out[idx] = in_1[idx] + in_2[idx];
}

//---------------------------------------------------------
void matrix_vector_product( const int m, const int n, const double matrix[m][n],
                            const double in[n], double out[m] )
{
  int idx, jdx;
  for( idx=0; idx<m; idx++ )
  {
    out[idx] = 0.0;
    for( jdx=0; jdx<n; jdx++ )
      out[idx] += matrix[idx][jdx] * in[jdx];
  }
}

//---------------------------------------------------------
void matrix_product( const int m, const int k, const int n,
                     const double A[m][k], const double B[k][n],
                     double C[m][n] )
{
  // A: m x k,  B: k x n  =>  C: m x n
  int idx, jdx, kdx;
  for(idx=0; idx<m; idx++)
  {
    for(jdx=0; jdx<n; jdx++)
    {
      C[idx][jdx] = 0.0;
      for(kdx = 0; kdx<k; kdx++)
        C[idx][jdx] += A[idx][kdx] * B[kdx][jdx];
    }
  }
}

//---------------------------------------------------------
double matrix_3by3_det( const double mat[3][3] )
{
  double det =
      mat[0][0] * ( mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1] )
    - mat[0][1] * ( mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0] )
    + mat[0][2] * ( mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0] );

  return det;
}

//---------------------------------------------------------
void matrix_t( const int m, const int n, const double A[m][n],
               double B[n][m] )
{
  int idx, jdx;
  for( idx=0; idx<m; idx++ )
  {
    for( jdx=0; jdx<n; jdx++ )
      B[jdx][idx] = A[idx][jdx];
  }
}

//---------------------------------------------------------
void skew_sym_3by3( const double in[3], double skew_sym[3][3] )
{
  int idx;
  for( idx=0; idx<3; idx++ )
    skew_sym[idx][idx] = 0.0;

  skew_sym[0][1] = -in[2];
  skew_sym[0][2] = in[1];
  skew_sym[1][2] = -in[0];

  skew_sym[1][0] = in[2];
  skew_sym[2][0] = -in[1];
  skew_sym[2][1] = in[0];
}

//---------------------------------------------------------
void mat_2_euler_angles( const double mat[3][3], double euler[3] )
{
  euler[0] = atan2( mat[2][1] , mat[2][2] );
  euler[1] = asin( -mat[2][0] );
  euler[2] = atan2( mat[1][0] , mat[0][0] );
}

//---------------------------------------------------------
void rotation_zyx( const double rot_z, const double rot_y,
                   const double rot_x, double R_zyx[3][3] )
{
  R_zyx[0][0] = cos( rot_z ) * cos( rot_y );
  R_zyx[0][1] = cos( rot_z ) * sin( rot_y ) * sin( rot_x ) -
    sin( rot_z ) * cos( rot_x );
  R_zyx[0][2] = cos( rot_z ) * sin( rot_y ) * cos( rot_x ) +
    sin( rot_z ) * sin( rot_x );
  R_zyx[1][0] = sin( rot_z ) * cos( rot_y );
  R_zyx[1][1] = sin( rot_z ) * sin( rot_y ) * sin( rot_x ) +
    cos( rot_z ) * cos( rot_x );
  R_zyx[1][2] = sin( rot_z ) * sin( rot_y ) * cos( rot_x ) -
    cos( rot_z ) * sin( rot_x );
  R_zyx[2][0] = -sin( rot_y );
  R_zyx[2][1] = cos( rot_y ) * sin( rot_x );
  R_zyx[2][2] = cos( rot_y ) * cos( rot_x );
}

//---------------------------------------------------------
void matrix_2_quaternions( const double mat[3][3], double q[4] )
{
  q[0] = 0.5 * sqrt( 1.0 + mat[0][0] + mat[1][1] + mat[2][2] );
  q[1] = 0.25 * ( mat[2][1] - mat[1][2] ) / q[0];
  q[2] = 0.25 * ( mat[0][2] - mat[2][0] ) / q[0];
  q[3] = 0.25 * ( mat[1][0] - mat[0][1] ) / q[0];
}

//---------------------------------------------------------
void quaternions_2_matrix( const double q[4], double mat[3][3] )
{
  mat[0][0] = 1.0 - 2.0 * ( q[2]* q[2] + q[3] * q[3] );
  mat[0][1] = 2.0 * ( q[1] * q[2] - q[0] * q[3] );
  mat[0][2] = 2.0 * ( q[1] * q[3] + q[0] * q[2] );

  mat[1][0] = 2.0 * ( q[1] * q[2] + q[0] * q[3] );
  mat[1][1] = 1.0 - 2.0 * ( q[1] * q[1] + q[3] * q[3] );
  mat[1][2] = 2.0 * ( q[2] * q[3] - q[0] * q[1] );

  mat[2][0] = 2.0 * ( q[1] * q[3] - q[0] * q[2] );
  mat[2][1] = 2.0 * ( q[2] * q[3] + q[0] * q[1] );
  mat[2][2] = 1.0 - 2.0 * ( q[1] * q[1] + q[2] * q[2] );
}

//---------------------------------------------------------
int invert3x3( double mat[3][3] )
{
  double inv[3][3];
  // Determinant via Laplace expansion
  double det =
    mat[0][0] * ( mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2] ) -
    mat[0][1] * ( mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0] ) +
    mat[0][2] * ( mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0] );

  int idx, jdx;
  double s = 0.0;
  for( idx=0; idx<3; idx++ )
  {
    for( jdx=0; jdx<3; jdx++ )
    {
      double a = fabs( mat[idx][jdx] );
      if( a > s )
        s = a;
    }
  }
  if( s == 0.0 )
    return 0;

  const double eps = 1e-12;
  if( fabs( det ) < eps * s * s * s )
    return 0;

  double invdet = 1.0 / det;

  // Calculate the adjugate matrix scaled by 1/det
  inv[0][0] = ( mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2] ) * invdet;
  inv[0][1] = ( mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2] ) * invdet;
  inv[0][2] = ( mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1] ) * invdet;

  inv[1][0] = ( mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2] ) * invdet;
  inv[1][1] = ( mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0] ) * invdet;
  inv[1][2] = ( mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2] ) * invdet;

  inv[2][0] = ( mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1] ) * invdet;
  inv[2][1] = ( mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1] ) * invdet;
  inv[2][2] = ( mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1] ) * invdet;

  memcpy( &mat[0][0], &inv[0][0], ( sizeof( double ) * 3 * 3 ) );

  return 1;
}

//---------------------------------------------------------
static int ludcmp( const int n, double A[n][n], int indx[n] )
{
  int i, j, k;
  int imax = 0;
  double vv[n];
  double d = 1.0;
  for( i=0; i<n; i++ )
  {
    double big = 0.0;
    for( j=0; j<n; j++ )
    {
      double tmp = fabs( A[j][i] );
      if( tmp > big )
        big = tmp;
    }

    if( big > 0.0 )
      vv[i] = 1.0 / big;
    else
      return 0;
  }

  for( j=0; j<n; j++ )
  {
    for( i=0; i<j; i++ )
    {
      double s = A[j][i];
      for( k=0; k<i; k++ )
        s -= A[k][i] * A[j][k];

      A[j][i] = s;
    }
    double big = 0.0;
    for( i=j; i<n; i++ )
    {
      double s = A[j][i];
      for( k=0; k<j; k++ )
        s -= A[k][i] * A[j][k];

      A[j][i] = s;
      double tmp = vv[i] * fabs( s );
      if( tmp >= big )
      {
        big = tmp;
        imax = i;
      }
    }
    if( j != imax )
    {
      for( k=0; k<n; k++ )
      {
        double tmp = A[k][imax];
        A[k][imax] = A[k][j];
        A[k][j] = tmp;
      }
      d = -d;
      vv[imax] = vv[j];
    }

    indx[j] = imax;
    if( A[j][j] == 0.0 )
      return 0;

    if( j != n-1 )
    {
      double tmp = 1.0 / A[j][j];
      for( i=j+1; i<n; i++ )
        A[j][i] *= tmp;
    }
  }

  return 1;
}

//---------------------------------------------------------
static void lubksb( const int n, const double A[n][n], const int indx[n], double b[n] )
{
  double s = 0;
  int idx, jdx;
  int ii = -1;

  for( idx=0; idx<n; idx++ )
  {
    int ip = indx[idx];
    s = b[ip];
    b[ip] = b[idx];

    if( ii>=0 )
    {
      for( jdx=ii; jdx<idx; jdx++ )
        s -= A[jdx][idx] * b[jdx];
    }
    else if ( s )
      ii = idx;

    b[idx] = s;
  }

  for( idx=n-1; idx>=0; idx-- )
  {
    s = b[idx];
    for( jdx=idx+1; jdx<n; jdx++ )
      s -= A[jdx][idx] * b[jdx];

    b[idx] = s / A[idx][idx];
  }
}

//---------------------------------------------------------
int lu_matrix_inv( const int n, double A[n][n] )
{
  double B[n][n];
  int indx[n];
  memcpy( B, A, ( sizeof( double ) * n * n ) );
  if( !ludcmp( n, B, indx ) )
    return 0;

  int idx, jdx;
  for( jdx=0; jdx<n; jdx++ )
  {
    for( idx=0; idx<n; idx++ )
      A[jdx][idx]=0.0;

    A[jdx][jdx]=1.0;
    lubksb( n, B, indx, A[jdx] );
  }

  return 1;
}
