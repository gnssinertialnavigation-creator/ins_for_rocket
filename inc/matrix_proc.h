#ifndef MATRIX_PROC_H_
#define MATRIX_PROC_H_

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "basic_const.h"

#ifdef __cplusplus
extern "C" {
#endif

void matrix_show( const int m, const int n, const double A[m][n] );
double vector_abs( const int n, const double in[n] );
void vector_norm( const int n, const double in[n], double out[n] );
void vector_sum( const int n, const double in_1[n], const double in_2[n], double out[n] );
void matrix_vector_product( const int m, const int n, const double matrix[m][n],
                            const double in[n], double out[m] );
void matrix_product( const int m, const int k, const int n,
                     const double A[m][k], const double B[k][n],
                     double C[m][n] );
void matrix_t( const int m, const int n, const double A[m][n],
               double B[n][m] );
void skew_sym_3by3( const double in[3], double skew_sym[3][3] );
double matrix_3by3_det( const double mat[3][3] );
void mat_2_euler_angles( const double mat[3][3], double euler[3] );
void rotation_zyx( const double rot_z, const double rot_y,
                   const double rot_x, double R_zyx[3][3] );
void matrix_2_quaternions( const double mat[3][3], double q[4] );
void quaternions_2_matrix( const double q[4], double mat[3][3] );
int invert3x3( double mat[3][3] );
int lu_matrix_inv( const int n, double A[n][n] );

#ifdef __cplusplus
}
#endif

#endif // #ifndef MATRIX_PROC_H_
