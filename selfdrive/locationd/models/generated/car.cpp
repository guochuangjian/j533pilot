#include "car.h"

namespace {
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 5.991464547107981;

/******************************************************************************
 *                      Code generated with sympy 1.7.1                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5826554309974982227) {
   out_5826554309974982227[0] = delta_x[0] + nom_x[0];
   out_5826554309974982227[1] = delta_x[1] + nom_x[1];
   out_5826554309974982227[2] = delta_x[2] + nom_x[2];
   out_5826554309974982227[3] = delta_x[3] + nom_x[3];
   out_5826554309974982227[4] = delta_x[4] + nom_x[4];
   out_5826554309974982227[5] = delta_x[5] + nom_x[5];
   out_5826554309974982227[6] = delta_x[6] + nom_x[6];
   out_5826554309974982227[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_1829394537637158322) {
   out_1829394537637158322[0] = -nom_x[0] + true_x[0];
   out_1829394537637158322[1] = -nom_x[1] + true_x[1];
   out_1829394537637158322[2] = -nom_x[2] + true_x[2];
   out_1829394537637158322[3] = -nom_x[3] + true_x[3];
   out_1829394537637158322[4] = -nom_x[4] + true_x[4];
   out_1829394537637158322[5] = -nom_x[5] + true_x[5];
   out_1829394537637158322[6] = -nom_x[6] + true_x[6];
   out_1829394537637158322[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_7631383496007094190) {
   out_7631383496007094190[0] = 1.0;
   out_7631383496007094190[1] = 0.0;
   out_7631383496007094190[2] = 0.0;
   out_7631383496007094190[3] = 0.0;
   out_7631383496007094190[4] = 0.0;
   out_7631383496007094190[5] = 0.0;
   out_7631383496007094190[6] = 0.0;
   out_7631383496007094190[7] = 0.0;
   out_7631383496007094190[8] = 0.0;
   out_7631383496007094190[9] = 1.0;
   out_7631383496007094190[10] = 0.0;
   out_7631383496007094190[11] = 0.0;
   out_7631383496007094190[12] = 0.0;
   out_7631383496007094190[13] = 0.0;
   out_7631383496007094190[14] = 0.0;
   out_7631383496007094190[15] = 0.0;
   out_7631383496007094190[16] = 0.0;
   out_7631383496007094190[17] = 0.0;
   out_7631383496007094190[18] = 1.0;
   out_7631383496007094190[19] = 0.0;
   out_7631383496007094190[20] = 0.0;
   out_7631383496007094190[21] = 0.0;
   out_7631383496007094190[22] = 0.0;
   out_7631383496007094190[23] = 0.0;
   out_7631383496007094190[24] = 0.0;
   out_7631383496007094190[25] = 0.0;
   out_7631383496007094190[26] = 0.0;
   out_7631383496007094190[27] = 1.0;
   out_7631383496007094190[28] = 0.0;
   out_7631383496007094190[29] = 0.0;
   out_7631383496007094190[30] = 0.0;
   out_7631383496007094190[31] = 0.0;
   out_7631383496007094190[32] = 0.0;
   out_7631383496007094190[33] = 0.0;
   out_7631383496007094190[34] = 0.0;
   out_7631383496007094190[35] = 0.0;
   out_7631383496007094190[36] = 1.0;
   out_7631383496007094190[37] = 0.0;
   out_7631383496007094190[38] = 0.0;
   out_7631383496007094190[39] = 0.0;
   out_7631383496007094190[40] = 0.0;
   out_7631383496007094190[41] = 0.0;
   out_7631383496007094190[42] = 0.0;
   out_7631383496007094190[43] = 0.0;
   out_7631383496007094190[44] = 0.0;
   out_7631383496007094190[45] = 1.0;
   out_7631383496007094190[46] = 0.0;
   out_7631383496007094190[47] = 0.0;
   out_7631383496007094190[48] = 0.0;
   out_7631383496007094190[49] = 0.0;
   out_7631383496007094190[50] = 0.0;
   out_7631383496007094190[51] = 0.0;
   out_7631383496007094190[52] = 0.0;
   out_7631383496007094190[53] = 0.0;
   out_7631383496007094190[54] = 1.0;
   out_7631383496007094190[55] = 0.0;
   out_7631383496007094190[56] = 0.0;
   out_7631383496007094190[57] = 0.0;
   out_7631383496007094190[58] = 0.0;
   out_7631383496007094190[59] = 0.0;
   out_7631383496007094190[60] = 0.0;
   out_7631383496007094190[61] = 0.0;
   out_7631383496007094190[62] = 0.0;
   out_7631383496007094190[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6042263353151296703) {
   out_6042263353151296703[0] = state[0];
   out_6042263353151296703[1] = state[1];
   out_6042263353151296703[2] = state[2];
   out_6042263353151296703[3] = state[3];
   out_6042263353151296703[4] = state[4];
   out_6042263353151296703[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6042263353151296703[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6042263353151296703[7] = state[7];
}
void F_fun(double *state, double dt, double *out_770756480427470884) {
   out_770756480427470884[0] = 1;
   out_770756480427470884[1] = 0;
   out_770756480427470884[2] = 0;
   out_770756480427470884[3] = 0;
   out_770756480427470884[4] = 0;
   out_770756480427470884[5] = 0;
   out_770756480427470884[6] = 0;
   out_770756480427470884[7] = 0;
   out_770756480427470884[8] = 0;
   out_770756480427470884[9] = 1;
   out_770756480427470884[10] = 0;
   out_770756480427470884[11] = 0;
   out_770756480427470884[12] = 0;
   out_770756480427470884[13] = 0;
   out_770756480427470884[14] = 0;
   out_770756480427470884[15] = 0;
   out_770756480427470884[16] = 0;
   out_770756480427470884[17] = 0;
   out_770756480427470884[18] = 1;
   out_770756480427470884[19] = 0;
   out_770756480427470884[20] = 0;
   out_770756480427470884[21] = 0;
   out_770756480427470884[22] = 0;
   out_770756480427470884[23] = 0;
   out_770756480427470884[24] = 0;
   out_770756480427470884[25] = 0;
   out_770756480427470884[26] = 0;
   out_770756480427470884[27] = 1;
   out_770756480427470884[28] = 0;
   out_770756480427470884[29] = 0;
   out_770756480427470884[30] = 0;
   out_770756480427470884[31] = 0;
   out_770756480427470884[32] = 0;
   out_770756480427470884[33] = 0;
   out_770756480427470884[34] = 0;
   out_770756480427470884[35] = 0;
   out_770756480427470884[36] = 1;
   out_770756480427470884[37] = 0;
   out_770756480427470884[38] = 0;
   out_770756480427470884[39] = 0;
   out_770756480427470884[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_770756480427470884[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_770756480427470884[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_770756480427470884[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_770756480427470884[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_770756480427470884[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_770756480427470884[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_770756480427470884[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_770756480427470884[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_770756480427470884[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_770756480427470884[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_770756480427470884[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_770756480427470884[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_770756480427470884[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_770756480427470884[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_770756480427470884[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_770756480427470884[56] = 0;
   out_770756480427470884[57] = 0;
   out_770756480427470884[58] = 0;
   out_770756480427470884[59] = 0;
   out_770756480427470884[60] = 0;
   out_770756480427470884[61] = 0;
   out_770756480427470884[62] = 0;
   out_770756480427470884[63] = 1;
}
void h_25(double *state, double *unused, double *out_5361102727716958993) {
   out_5361102727716958993[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3207636012344091342) {
   out_3207636012344091342[0] = 0;
   out_3207636012344091342[1] = 0;
   out_3207636012344091342[2] = 0;
   out_3207636012344091342[3] = 0;
   out_3207636012344091342[4] = 0;
   out_3207636012344091342[5] = 0;
   out_3207636012344091342[6] = 1;
   out_3207636012344091342[7] = 0;
}
void h_24(double *state, double *unused, double *out_2493276794756230534) {
   out_2493276794756230534[0] = state[4];
   out_2493276794756230534[1] = state[5];
}
void H_24(double *state, double *unused, double *out_3137807261716376900) {
   out_3137807261716376900[0] = 0;
   out_3137807261716376900[1] = 0;
   out_3137807261716376900[2] = 0;
   out_3137807261716376900[3] = 0;
   out_3137807261716376900[4] = 1;
   out_3137807261716376900[5] = 0;
   out_3137807261716376900[6] = 0;
   out_3137807261716376900[7] = 0;
   out_3137807261716376900[8] = 0;
   out_3137807261716376900[9] = 0;
   out_3137807261716376900[10] = 0;
   out_3137807261716376900[11] = 0;
   out_3137807261716376900[12] = 0;
   out_3137807261716376900[13] = 1;
   out_3137807261716376900[14] = 0;
   out_3137807261716376900[15] = 0;
}
void h_30(double *state, double *unused, double *out_4902901365558397331) {
   out_4902901365558397331[0] = state[4];
}
void H_30(double *state, double *unused, double *out_6406262898605091938) {
   out_6406262898605091938[0] = 0;
   out_6406262898605091938[1] = 0;
   out_6406262898605091938[2] = 0;
   out_6406262898605091938[3] = 0;
   out_6406262898605091938[4] = 1;
   out_6406262898605091938[5] = 0;
   out_6406262898605091938[6] = 0;
   out_6406262898605091938[7] = 0;
}
void h_26(double *state, double *unused, double *out_3893408098707599809) {
   out_3893408098707599809[0] = state[7];
}
void H_26(double *state, double *unused, double *out_8912229167283300615) {
   out_8912229167283300615[0] = 0;
   out_8912229167283300615[1] = 0;
   out_8912229167283300615[2] = 0;
   out_8912229167283300615[3] = 0;
   out_8912229167283300615[4] = 0;
   out_8912229167283300615[5] = 0;
   out_8912229167283300615[6] = 0;
   out_8912229167283300615[7] = 1;
}
void h_27(double *state, double *unused, double *out_5964418612393436510) {
   out_5964418612393436510[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7693844886441717250) {
   out_7693844886441717250[0] = 0;
   out_7693844886441717250[1] = 0;
   out_7693844886441717250[2] = 0;
   out_7693844886441717250[3] = 1;
   out_7693844886441717250[4] = 0;
   out_7693844886441717250[5] = 0;
   out_7693844886441717250[6] = 0;
   out_7693844886441717250[7] = 0;
}
void h_29(double *state, double *unused, double *out_1866115314114064646) {
   out_1866115314114064646[0] = state[1];
}
void H_29(double *state, double *unused, double *out_5831704284749555763) {
   out_5831704284749555763[0] = 0;
   out_5831704284749555763[1] = 1;
   out_5831704284749555763[2] = 0;
   out_5831704284749555763[3] = 0;
   out_5831704284749555763[4] = 0;
   out_5831704284749555763[5] = 0;
   out_5831704284749555763[6] = 0;
   out_5831704284749555763[7] = 0;
}
void h_28(double *state, double *unused, double *out_3369080615822404081) {
   out_3369080615822404081[0] = state[5];
   out_3369080615822404081[1] = state[6];
}
void H_28(double *state, double *unused, double *out_8291104333756679850) {
   out_8291104333756679850[0] = 0;
   out_8291104333756679850[1] = 0;
   out_8291104333756679850[2] = 0;
   out_8291104333756679850[3] = 0;
   out_8291104333756679850[4] = 0;
   out_8291104333756679850[5] = 1;
   out_8291104333756679850[6] = 0;
   out_8291104333756679850[7] = 0;
   out_8291104333756679850[8] = 0;
   out_8291104333756679850[9] = 0;
   out_8291104333756679850[10] = 0;
   out_8291104333756679850[11] = 0;
   out_8291104333756679850[12] = 0;
   out_8291104333756679850[13] = 0;
   out_8291104333756679850[14] = 1;
   out_8291104333756679850[15] = 0;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5826554309974982227) {
  err_fun(nom_x, delta_x, out_5826554309974982227);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1829394537637158322) {
  inv_err_fun(nom_x, true_x, out_1829394537637158322);
}
void car_H_mod_fun(double *state, double *out_7631383496007094190) {
  H_mod_fun(state, out_7631383496007094190);
}
void car_f_fun(double *state, double dt, double *out_6042263353151296703) {
  f_fun(state,  dt, out_6042263353151296703);
}
void car_F_fun(double *state, double dt, double *out_770756480427470884) {
  F_fun(state,  dt, out_770756480427470884);
}
void car_h_25(double *state, double *unused, double *out_5361102727716958993) {
  h_25(state, unused, out_5361102727716958993);
}
void car_H_25(double *state, double *unused, double *out_3207636012344091342) {
  H_25(state, unused, out_3207636012344091342);
}
void car_h_24(double *state, double *unused, double *out_2493276794756230534) {
  h_24(state, unused, out_2493276794756230534);
}
void car_H_24(double *state, double *unused, double *out_3137807261716376900) {
  H_24(state, unused, out_3137807261716376900);
}
void car_h_30(double *state, double *unused, double *out_4902901365558397331) {
  h_30(state, unused, out_4902901365558397331);
}
void car_H_30(double *state, double *unused, double *out_6406262898605091938) {
  H_30(state, unused, out_6406262898605091938);
}
void car_h_26(double *state, double *unused, double *out_3893408098707599809) {
  h_26(state, unused, out_3893408098707599809);
}
void car_H_26(double *state, double *unused, double *out_8912229167283300615) {
  H_26(state, unused, out_8912229167283300615);
}
void car_h_27(double *state, double *unused, double *out_5964418612393436510) {
  h_27(state, unused, out_5964418612393436510);
}
void car_H_27(double *state, double *unused, double *out_7693844886441717250) {
  H_27(state, unused, out_7693844886441717250);
}
void car_h_29(double *state, double *unused, double *out_1866115314114064646) {
  h_29(state, unused, out_1866115314114064646);
}
void car_H_29(double *state, double *unused, double *out_5831704284749555763) {
  H_29(state, unused, out_5831704284749555763);
}
void car_h_28(double *state, double *unused, double *out_3369080615822404081) {
  h_28(state, unused, out_3369080615822404081);
}
void car_H_28(double *state, double *unused, double *out_8291104333756679850) {
  H_28(state, unused, out_8291104333756679850);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
