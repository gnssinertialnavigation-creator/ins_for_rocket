#include "basic_const.h"
#include "matrix_proc.h"
#include "navigation_process_function.h"
#include "ins_01.h"
#include "ins_02.h"
#include "ins_03.h"
#include "ins_04.h"
#include "rocket_body_force_10hz.h"
#include "rocket_body_omega_10hz.h"

static void rocket_debuger( FILE *rocket_fp, double tsince, DateTime* utc_time_p,
                            double p_ecef[3], double v_ecef[3], double qu_b2e[4],
                            double eular_b2l[3] );
static void ins_debuger( FILE *ins_01_fp, FILE *ins_02_fp,
                         FILE *ins_03_fp, FILE *ins_04_fp,
                         const double tsince, const DateTime* utc_time_p,
                         const double p_ecef[3], const double v_ecef[3],
                         const double qu_b2e[4], const double eular_b2l[3] );
static void kml_head( FILE *fp );
static void kml_tail( FILE *fp );
static void lla_2_kml( FILE *fp, double lat, double lon, double alt );

//---------------------------------------------------------
static void rocket_debuger( FILE *rocket_fp, double tsince, DateTime* utc_time_p,
                            double p_ecef[3], double v_ecef[3], double qu_b2e[4],
                            double eular_b2l[3] )
{
  fprintf( rocket_fp, "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,"
           "ecef_v,%.10f,%.10f,%.10f,qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "eular_b2l,%.10f,%.10f,%.10f\n",
           tsince, utc_time_p->year, utc_time_p->month, utc_time_p->day, utc_time_p->hours,
           p_ecef[0] * 1E-3, p_ecef[1] * 1E-3, p_ecef[2] * 1E-3,
           v_ecef[0] * 1E-3, v_ecef[1] * 1E-3, v_ecef[2] * 1E-3,
           qu_b2e[0], qu_b2e[1], qu_b2e[2], qu_b2e[3], vector_abs( 4, qu_b2e ),
           eular_b2l[0] * RAD_2_DEG, eular_b2l[1] * RAD_2_DEG, eular_b2l[2] * RAD_2_DEG );
}

//---------------------------------------------------------
static void ins_debuger( FILE *ins_01_fp, FILE *ins_02_fp,
                         FILE *ins_03_fp, FILE *ins_04_fp,
                         const double tsince, const DateTime* utc_time_p,
                         const double p_ecef[3], const double v_ecef[3],
                         const double qu_b2e[4], const double eular_b2l[3] )
{
  double pos[3], vel[3], qu[4], bf[3], bg[3];
  get_ins_01_pos( &pos[0] );
  get_ins_01_vel( &vel[0] );
  get_ins_01_qu_b2e( &qu[0] );
  get_ins_01_bias_force( &bf[0] );
  get_ins_01_bias_gyro( &bg[0] );
  double ins_lat, ins_lon, ins_height;
  double ins_eular_b2l[3], C_e2l[3][3], C_b2e[3][3], C_b2l[3][3];
  ecef2lla( pos, &ins_lat, &ins_lon, &ins_height );
  ecef2enu_mat( ins_lat, ins_lon, C_e2l );
  quaternions_2_matrix( &qu[0], C_b2e );
  matrix_product( 3, 3, 3, C_e2l, C_b2e, C_b2l );
  mat_2_euler_angles( C_b2l, ins_eular_b2l );

  fprintf( ins_01_fp, "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "ecef_v,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "eular_b2l,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "bf,%.10f,%.10f,%.10f,bg,%.10f,%.10f,%.10f\n",
           tsince, utc_time_p->year, utc_time_p->month, utc_time_p->day, utc_time_p->hours,
           pos[0] * 1E-3, pos[1] * 1E-3, pos[2] * 1E-3,
           ( p_ecef[0] - pos[0] ) * 1E-3, ( p_ecef[1] - pos[1] ) * 1E-3,
           ( p_ecef[2] - pos[2] ) * 1E-3,
           vel[0] * 1E-3, vel[1] * 1E-3, vel[2] * 1E-3,
           ( v_ecef[0] - vel[0] ) * 1E-3, ( v_ecef[1] - vel[1] ) * 1E-3,
           ( v_ecef[2] - vel[2] ) * 1E-3,
           qu[0], qu[1], qu[2], qu[3], vector_abs( 4, &qu[0] ),
           qu_b2e[0] - qu[0], qu_b2e[1] - qu[1], qu_b2e[2] - qu[2], qu_b2e[3] - qu[3],
           vector_abs( 4, &qu_b2e[0] ),
           ins_eular_b2l[0] * RAD_2_DEG, ins_eular_b2l[1] * RAD_2_DEG, ins_eular_b2l[2] * RAD_2_DEG,
           ( eular_b2l[0] - ins_eular_b2l[0] ) * RAD_2_DEG,
           ( eular_b2l[1] - ins_eular_b2l[1] ) * RAD_2_DEG,
           ( eular_b2l[2] - ins_eular_b2l[2] ) * RAD_2_DEG,
           bf[0] * N_2_G, bf[1] * N_2_G, bf[2] * N_2_G,
           bg[0] * RAD_2_DEG, bg[1] * RAD_2_DEG, bg[2] * RAD_2_DEG );

  get_ins_02_pos( &pos[0] );
  get_ins_02_vel( &vel[0] );
  get_ins_02_qu_b2e( &qu[0] );
  get_ins_02_bias_force( &bf[0] );
  get_ins_02_bias_gyro( &bg[0] );
  ecef2lla( pos, &ins_lat, &ins_lon, &ins_height );
  ecef2enu_mat( ins_lat, ins_lon, C_e2l );
  quaternions_2_matrix( &qu[0], C_b2e );
  matrix_product( 3, 3, 3, C_e2l, C_b2e, C_b2l );
  mat_2_euler_angles( C_b2l, ins_eular_b2l );

  fprintf( ins_02_fp, "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "ecef_v,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "eular_b2l,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "bf,%.10f,%.10f,%.10f,bg,%.10f,%.10f,%.10f\n",
           tsince, utc_time_p->year, utc_time_p->month, utc_time_p->day, utc_time_p->hours,
           pos[0] * 1E-3, pos[1] * 1E-3, pos[2] * 1E-3,
           ( p_ecef[0] - pos[0] ) * 1E-3, ( p_ecef[1] - pos[1] ) * 1E-3,
           ( p_ecef[2] - pos[2] ) * 1E-3,
           vel[0] * 1E-3, vel[1] * 1E-3, vel[2] * 1E-3,
           ( v_ecef[0] - vel[0] ) * 1E-3, ( v_ecef[1] - vel[1] ) * 1E-3,
           ( v_ecef[2] - vel[2] ) * 1E-3,
           qu[0], qu[1], qu[2], qu[3], vector_abs( 4, &qu[0] ),
           qu_b2e[0] - qu[0], qu_b2e[1] - qu[1], qu_b2e[2] - qu[2], qu_b2e[3] - qu[3],
           vector_abs( 4, &qu_b2e[0] ),
           ins_eular_b2l[0] * RAD_2_DEG, ins_eular_b2l[1] * RAD_2_DEG, ins_eular_b2l[2] * RAD_2_DEG,
           ( eular_b2l[0] - ins_eular_b2l[0] ) * RAD_2_DEG,
           ( eular_b2l[1] - ins_eular_b2l[1] ) * RAD_2_DEG,
           ( eular_b2l[2] - ins_eular_b2l[2] ) * RAD_2_DEG,
           bf[0] * N_2_G, bf[1] * N_2_G, bf[2] * N_2_G,
           bg[0] * RAD_2_DEG, bg[1] * RAD_2_DEG, bg[2] * RAD_2_DEG );

  get_ins_03_pos( &pos[0] );
  get_ins_03_vel( &vel[0] );
  get_ins_03_qu_b2e( &qu[0] );
  get_ins_03_bias_force( &bf[0] );
  get_ins_03_bias_gyro( &bg[0] );
  ecef2lla( pos, &ins_lat, &ins_lon, &ins_height );
  ecef2enu_mat( ins_lat, ins_lon, C_e2l );
  quaternions_2_matrix( &qu[0], C_b2e );
  matrix_product( 3, 3, 3, C_e2l, C_b2e, C_b2l );
  mat_2_euler_angles( C_b2l, ins_eular_b2l );

  fprintf( ins_03_fp, "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "ecef_v,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "eular_b2l,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "bf,%.10f,%.10f,%.10f,bg,%.10f,%.10f,%.10f\n",
           tsince, utc_time_p->year, utc_time_p->month, utc_time_p->day, utc_time_p->hours,
           pos[0] * 1E-3, pos[1] * 1E-3, pos[2] * 1E-3,
           ( p_ecef[0] - pos[0] ) * 1E-3, ( p_ecef[1] - pos[1] ) * 1E-3,
           ( p_ecef[2] - pos[2] ) * 1E-3,
           vel[0] * 1E-3, vel[1] * 1E-3, vel[2] * 1E-3,
           ( v_ecef[0] - vel[0] ) * 1E-3, ( v_ecef[1] - vel[1] ) * 1E-3,
           ( v_ecef[2] - vel[2] ) * 1E-3,
           qu[0], qu[1], qu[2], qu[3], vector_abs( 4, &qu[0] ),
           qu_b2e[0] - qu[0], qu_b2e[1] - qu[1], qu_b2e[2] - qu[2], qu_b2e[3] - qu[3],
           vector_abs( 4, &qu_b2e[0] ),
           ins_eular_b2l[0] * RAD_2_DEG, ins_eular_b2l[1] * RAD_2_DEG, ins_eular_b2l[2] * RAD_2_DEG,
           ( eular_b2l[0] - ins_eular_b2l[0] ) * RAD_2_DEG,
           ( eular_b2l[1] - ins_eular_b2l[1] ) * RAD_2_DEG,
           ( eular_b2l[2] - ins_eular_b2l[2] ) * RAD_2_DEG,
           bf[0] * N_2_G, bf[1] * N_2_G, bf[2] * N_2_G,
           bg[0] * RAD_2_DEG, bg[1] * RAD_2_DEG, bg[2] * RAD_2_DEG );

  get_ins_04_pos( &pos[0] );
  get_ins_04_vel( &vel[0] );
  get_ins_04_qu_b2e( &qu[0] );
  get_ins_04_bias_force( &bf[0] );
  get_ins_04_bias_gyro( &bg[0] );
  ecef2lla( pos, &ins_lat, &ins_lon, &ins_height );
  ecef2enu_mat( ins_lat, ins_lon, C_e2l );
  quaternions_2_matrix( &qu[0], C_b2e );
  matrix_product( 3, 3, 3, C_e2l, C_b2e, C_b2l );
  mat_2_euler_angles( C_b2l, ins_eular_b2l );

  fprintf( ins_04_fp, "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "ecef_v,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "eular_b2l,%.10f,%.10f,%.10f,%.10f,%.10f,%.10f,"
           "bf,%.10f,%.10f,%.10f,bg,%.10f,%.10f,%.10f\n",
           tsince, utc_time_p->year, utc_time_p->month, utc_time_p->day, utc_time_p->hours,
           pos[0] * 1E-3, pos[1] * 1E-3, pos[2] * 1E-3,
           ( p_ecef[0] - pos[0] ) * 1E-3, ( p_ecef[1] - pos[1] ) * 1E-3,
           ( p_ecef[2] - pos[2] ) * 1E-3,
           vel[0] * 1E-3, vel[1] * 1E-3, vel[2] * 1E-3,
           ( v_ecef[0] - vel[0] ) * 1E-3, ( v_ecef[1] - vel[1] ) * 1E-3,
           ( v_ecef[2] - vel[2] ) * 1E-3,
           qu[0], qu[1], qu[2], qu[3], vector_abs( 4, &qu[0] ),
           qu_b2e[0] - qu[0], qu_b2e[1] - qu[1], qu_b2e[2] - qu[2], qu_b2e[3] - qu[3],
           vector_abs( 4, &qu_b2e[0] ),
           ins_eular_b2l[0] * RAD_2_DEG, ins_eular_b2l[1] * RAD_2_DEG, ins_eular_b2l[2] * RAD_2_DEG,
           ( eular_b2l[0] - ins_eular_b2l[0] ) * RAD_2_DEG,
           ( eular_b2l[1] - ins_eular_b2l[1] ) * RAD_2_DEG,
           ( eular_b2l[2] - ins_eular_b2l[2] ) * RAD_2_DEG,
           bf[0] * N_2_G, bf[1] * N_2_G, bf[2] * N_2_G,
           bg[0] * RAD_2_DEG, bg[1] * RAD_2_DEG, bg[2] * RAD_2_DEG );
}

//---------------------------------------------------------
static void kml_head( FILE *fp )
{
  fprintf( fp, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" );
  fprintf( fp, "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n" );
  fprintf( fp, "<Document>\n" );
  fprintf( fp, "  <Style id=\"greenDot\">\n" );
  fprintf( fp, "    <IconStyle><color>ff00ff00</color></IconStyle>\n" );
  fprintf( fp, "  </Style>\n" );
}

//---------------------------------------------------------
static void kml_tail( FILE *fp )
{
  fprintf( fp, "</Document>\n" );
  fprintf( fp, "</kml>\n" );
}

//---------------------------------------------------------
static void lla_2_kml( FILE *fp, double lat, double lon, double alt )
{
  fprintf( fp, "  <Placemark>\n" );
  fprintf( fp, "    <styleUrl>#greenDot</styleUrl>\n" );
  fprintf( fp, "    <Point>\n" );
  fprintf( fp, "      <altitudeMode>absolute</altitudeMode>\n" );
  fprintf( fp, "      <coordinates>%.10f,%.10f,%.10f</coordinates>\n",
           lon, lat, alt );
  fprintf( fp, "    </Point>\n" );
  fprintf( fp, "  </Placemark>\n" );
}

//---------------------------------------------------------
int main()
{
  FILE *rocket_fp = fopen("rocket.txt", "w");
  FILE *ins_01_fp = fopen("ins_01.txt", "w");
  FILE *ins_02_fp = fopen("ins_02.txt", "w");
  FILE *ins_03_fp = fopen("ins_03.txt", "w");
  FILE *ins_04_fp = fopen("ins_04.txt", "w");
  FILE *kml_fp = fopen("rocket.kml", "w");
  // Write KML header
  kml_head( kml_fp );

  double tsince = -900.0;
  const double tsince_end = 900.0;
  double p_ecef[3] = { 0, 0, 0 };
  double v_ecef[3] = { 0, 0, 0 };
  double p_eci[3] = { 0, 0, 0 };
  double v_eci[3] = { 0, 0, 0 };
  double lat = 30.4013 * DEG_2_RAD;
  double lon = 130.9702 * DEG_2_RAD;
  double height = 18.0;
  lla2ecef( lat, lon, height, p_ecef );
  double tilt_x = -6 * DEG_2_RAD;
  double tilt_y = 8 * DEG_2_RAD;
  double lauch_azimuth = 85 * DEG_2_RAD; // pointing of body x.
  double C_b2l[3][3];
  rotation_zyx( lauch_azimuth, tilt_y, tilt_x, C_b2l );
  double eular_b2l[3];
  double v_enu[3];
  double C_b2e[3][3];
  double C_e2l[3][3];
  double C_l2e[3][3];
  double abs_v = 0;
  ecef2enu_mat( lat, lon, C_e2l );
  matrix_t( 3, 3, C_e2l, C_l2e );
  matrix_product( 3, 3, 3, C_l2e, C_b2l, C_b2e );
  double qu_b2e[4];
  matrix_2_quaternions( C_b2e, qu_b2e );

  double f0_body[3] = { 0, 0, 0 };
  double gyro0_body[3] = { 0, 0, 0 };
  double f_body[3] = { 0, 0, 0 };
  double gyro_body[3] = { 0, 0, 0 };

  DateTime utc_time;
  utc_time.year = 2027;
  utc_time.month = 12;
  utc_time.day = 31;
  utc_time.hours = 23.95;
  ecef_to_eci( &utc_time, p_ecef, v_ecef, p_eci, v_eci );

  int data_hz = 10;
  double dt = 1.0 / data_hz;

  double p_ecef_error[3];
  double v_ecef_error[3];
  double qu_b2e_error[4];
  p_ecef_error[0] = p_ecef[0] + 100.2;
  p_ecef_error[1] = p_ecef[1] + 50.1;
  p_ecef_error[2] = p_ecef[2] - 80.8;
  v_ecef_error[0] = v_ecef[0] + 0.2;
  v_ecef_error[1] = v_ecef[1] - 0.3;
  v_ecef_error[2] = v_ecef[2] + 0.9;
  double C_b2e_error[3][3], C_skew[3][3], C_b2l_error[3][3];
  //rotation_zyx( 0 * DEG_2_RAD, 0 * DEG_2_RAD, 0 * DEG_2_RAD, C_skew );
  rotation_zyx( 26 * DEG_2_RAD, -2.1 * DEG_2_RAD, 4.5 * DEG_2_RAD, C_skew );
  matrix_product( 3, 3, 3, C_skew, C_b2l, C_b2l_error );
  matrix_product( 3, 3, 3, C_l2e, C_b2l_error, C_b2e_error );
  matrix_2_quaternions( C_b2e_error, qu_b2e_error );
  reset_ins_01( p_ecef_error, v_ecef_error, qu_b2e_error );
  reset_ins_02( p_ecef_error, v_ecef_error, qu_b2e_error );
  reset_ins_03( p_ecef_error, v_ecef_error, qu_b2e_error );
  reset_ins_04( p_ecef_error, v_ecef_error, qu_b2e_error );
  double bf_01[3] = { -2.3 * imu_force_bias_std_01, 4.3 * imu_force_bias_std_01, 5.4 * imu_force_bias_std_01 };
  double bg_01[3] = { -3.2 * imu_gyro_bias_std_01, -3.4 * imu_gyro_bias_std_01, 4.5 * imu_gyro_bias_std_01 };
  double bf_02[3] = { -2.3 * imu_force_bias_std_02, 4.3 * imu_force_bias_std_02, 5.4 * imu_force_bias_std_02 };
  double bg_02[3] = { -3.2 * imu_gyro_bias_std_02, -3.4 * imu_gyro_bias_std_02, 4.5 * imu_gyro_bias_std_02 };
  double bf_03[3] = { -2.3 * imu_force_bias_std_03, 4.3 * imu_force_bias_std_03, 5.4 * imu_force_bias_std_03 };
  double bg_03[3] = { -3.2 * imu_gyro_bias_std_03, -3.4 * imu_gyro_bias_std_03, 4.5 * imu_gyro_bias_std_03 };
  double bf_04[3] = { -2.3 * imu_force_bias_std_04, 4.3 * imu_force_bias_std_04, 5.4 * imu_force_bias_std_04 };
  double bg_04[3] = { -3.2 * imu_gyro_bias_std_04, -3.4 * imu_gyro_bias_std_04, 4.5 * imu_gyro_bias_std_04 };

  int idx, jdx;
  idx = 0;
  while( tsince < tsince_end )
  {
    double time_last = tsince;
    double time_now = time_last + dt;

    get_rocket_body_force( time_last, C_b2e, p_ecef, f0_body );
    get_rocket_body_force( time_now, C_b2e, p_ecef, f_body );

    get_rocket_body_omega( time_last, C_b2e, gyro0_body );
    get_rocket_body_omega( time_now, C_b2e, gyro_body );

    ecef2lla( p_ecef, &lat, &lon, &height );
    ecef2enu( lat, lon, C_e2l, v_ecef, v_enu );
    matrix_product( 3, 3, 3, C_e2l, C_b2e, C_b2l );
    mat_2_euler_angles( C_b2l, eular_b2l );
    abs_v = v_enu[0] * v_enu[0] + v_enu[1] * v_enu[1] + v_enu[2] * v_enu[2];
    abs_v = sqrt( abs_v );
    if( height < 0 )
      break;

    // 1 sec printf.
    if( ( idx % data_hz ) == 0 )
    {
      lla_2_kml( kml_fp, lat * RAD_2_DEG, lon * RAD_2_DEG, height );
      rocket_debuger( rocket_fp, time_last, &utc_time, p_ecef, v_ecef, qu_b2e, eular_b2l );
      ins_debuger( ins_01_fp, ins_02_fp, ins_03_fp, ins_04_fp,
                   time_last, &utc_time, p_ecef, v_ecef, qu_b2e, eular_b2l );
      if( DEBUGER_ON )
      {
        printf( "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,ecef_v,%.10f,%.10f,%.10f,"
                "lat,%.10f,lon,%.10f,h,%.10f,enu_v,%.10f,%.10f,%.10f,abs_v,%.10f,"
                "qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,eular_b2l,%.10f,%.10f,%.10f,"
                "det_C_b2e,%.10f\n",
                time_last, utc_time.year, utc_time.month, utc_time.day, utc_time.hours,
                p_ecef[0] * 1E-3, p_ecef[1] * 1E-3, p_ecef[2] * 1E-3,
                v_ecef[0] * 1E-3, v_ecef[1] * 1E-3, v_ecef[2] * 1E-3,
                lat * RAD_2_DEG, lon * RAD_2_DEG, height,
                v_enu[0] * 1E-3, v_enu[1] * 1E-3, v_enu[2] * 1E-3, abs_v * 1E-3,
                qu_b2e[0], qu_b2e[1], qu_b2e[2], qu_b2e[3], vector_abs( 4, qu_b2e ),
                eular_b2l[0] * RAD_2_DEG, eular_b2l[1] * RAD_2_DEG, eular_b2l[2] * RAD_2_DEG,
                matrix_3by3_det( C_b2e ) );
      }
    }
    // true rocket states
    rk4_states_propagator( time_last, p_ecef, v_ecef, qu_b2e,
                           f0_body, gyro0_body, f_body, gyro_body, dt );
    quaternions_2_matrix( &qu_b2e[0], C_b2e );
    add_seconds_to_utc( &utc_time, dt );
    ecef_to_eci( &utc_time, p_ecef, v_ecef, p_eci, v_eci );

    // ins
    double f0_b[3], gy0_b[3];
    for( jdx=0; jdx<3; jdx++ )
    {
      f0_b[jdx] = f0_body[jdx] + bf_01[jdx];
      gy0_b[jdx] = gyro0_body[jdx] + bg_01[jdx];
    }
    ins_01( time_last, p_ecef, v_ecef, C_b2e, f0_b, gy0_b, dt );
    for( jdx=0; jdx<3; jdx++ )
    {
      f0_b[jdx] = f0_body[jdx] + bf_02[jdx];
      gy0_b[jdx] = gyro0_body[jdx] + bg_02[jdx];
    }
    ins_02( time_last, p_ecef, v_ecef, C_b2e, f0_b, gy0_b, dt );
    for( jdx=0; jdx<3; jdx++ )
    {
      f0_b[jdx] = f0_body[jdx] + bf_03[jdx];
      gy0_b[jdx] = gyro0_body[jdx] + bg_03[jdx];
    }
    ins_03( time_last, p_ecef, v_ecef, qu_b2e, f0_b, gy0_b, dt );
    for( jdx=0; jdx<3; jdx++ )
    {
      f0_b[jdx] = f0_body[jdx] + bf_04[jdx];
      gy0_b[jdx] = gyro0_body[jdx] + bg_04[jdx];
    }
    ins_04( time_last, p_ecef, v_ecef, qu_b2e, f0_b, gy0_b, dt );

    tsince = time_now;
    idx++;
  }

  // last printf
  ecef2lla( p_ecef, &lat, &lon, &height );
  ecef2enu( lat, lon, C_e2l, v_ecef, v_enu );
  matrix_product( 3, 3, 3, C_e2l, C_b2e, C_b2l );
  mat_2_euler_angles( C_b2l, eular_b2l );
  abs_v = v_enu[0] * v_enu[0] + v_enu[1] * v_enu[1] + v_enu[2] * v_enu[2];
  abs_v = sqrt( abs_v );
  lla_2_kml( kml_fp, lat * RAD_2_DEG, lon * RAD_2_DEG, height );
  rocket_debuger( rocket_fp, tsince, &utc_time, p_ecef, v_ecef, qu_b2e, eular_b2l );
  ins_debuger( ins_01_fp, ins_02_fp, ins_03_fp, ins_04_fp,
               tsince, &utc_time, p_ecef, v_ecef, qu_b2e, eular_b2l );
  if( DEBUGER_ON )
  {
    printf( "tsince,%lf,%d,%d,%d,%lf,ecef_p,%.10f,%.10f,%.10f,ecef_v,%.10f,%.10f,%.10f,"
            "lat,%.10f,lon,%.10f,h,%.10f,enu_v,%.10f,%.10f,%.10f,abs_v,%.10f,"
            "qu_b2e,%.10f,%.10f,%.10f,%.10f,%.10f,eular_b2l,%.10f,%.10f,%.10f,"
            "det_C_b2e,%.10f\n",
            tsince, utc_time.year, utc_time.month, utc_time.day, utc_time.hours,
            p_ecef[0] * 1E-3, p_ecef[1] * 1E-3, p_ecef[2] * 1E-3,
            v_ecef[0] * 1E-3, v_ecef[1] * 1E-3, v_ecef[2] * 1E-3,
            lat * RAD_2_DEG, lon * RAD_2_DEG, height,
            v_enu[0] * 1E-3, v_enu[1] * 1E-3, v_enu[2] * 1E-3, abs_v * 1E-3,
            qu_b2e[0], qu_b2e[1], qu_b2e[2], qu_b2e[3], vector_abs( 4, qu_b2e ),
            eular_b2l[0] * RAD_2_DEG, eular_b2l[1] * RAD_2_DEG, eular_b2l[2] * RAD_2_DEG,
            matrix_3by3_det( C_b2e ) );
  }

  fclose( rocket_fp );
  fclose( ins_01_fp );
  fclose( ins_02_fp );
  fclose( ins_03_fp );
  fclose( ins_04_fp );
  kml_tail( kml_fp );
  fclose( kml_fp );

  return 0;
}
