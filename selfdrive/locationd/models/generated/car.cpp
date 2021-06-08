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
void err_fun(double *nom_x, double *delta_x, double *out_3837833141682255610) {
   out_3837833141682255610[0] = delta_x[0] + nom_x[0];
   out_3837833141682255610[1] = delta_x[1] + nom_x[1];
   out_3837833141682255610[2] = delta_x[2] + nom_x[2];
   out_3837833141682255610[3] = delta_x[3] + nom_x[3];
   out_3837833141682255610[4] = delta_x[4] + nom_x[4];
   out_3837833141682255610[5] = delta_x[5] + nom_x[5];
   out_3837833141682255610[6] = delta_x[6] + nom_x[6];
   out_3837833141682255610[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_3068960279000246371) {
   out_3068960279000246371[0] = -nom_x[0] + true_x[0];
   out_3068960279000246371[1] = -nom_x[1] + true_x[1];
   out_3068960279000246371[2] = -nom_x[2] + true_x[2];
   out_3068960279000246371[3] = -nom_x[3] + true_x[3];
   out_3068960279000246371[4] = -nom_x[4] + true_x[4];
   out_3068960279000246371[5] = -nom_x[5] + true_x[5];
   out_3068960279000246371[6] = -nom_x[6] + true_x[6];
   out_3068960279000246371[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_3240948688174513378) {
   out_3240948688174513378[0] = 1.0;
   out_3240948688174513378[1] = 0.0;
   out_3240948688174513378[2] = 0.0;
   out_3240948688174513378[3] = 0.0;
   out_3240948688174513378[4] = 0.0;
   out_3240948688174513378[5] = 0.0;
   out_3240948688174513378[6] = 0.0;
   out_3240948688174513378[7] = 0.0;
   out_3240948688174513378[8] = 0.0;
   out_3240948688174513378[9] = 1.0;
   out_3240948688174513378[10] = 0.0;
   out_3240948688174513378[11] = 0.0;
   out_3240948688174513378[12] = 0.0;
   out_3240948688174513378[13] = 0.0;
   out_3240948688174513378[14] = 0.0;
   out_3240948688174513378[15] = 0.0;
   out_3240948688174513378[16] = 0.0;
   out_3240948688174513378[17] = 0.0;
   out_3240948688174513378[18] = 1.0;
   out_3240948688174513378[19] = 0.0;
   out_3240948688174513378[20] = 0.0;
   out_3240948688174513378[21] = 0.0;
   out_3240948688174513378[22] = 0.0;
   out_3240948688174513378[23] = 0.0;
   out_3240948688174513378[24] = 0.0;
   out_3240948688174513378[25] = 0.0;
   out_3240948688174513378[26] = 0.0;
   out_3240948688174513378[27] = 1.0;
   out_3240948688174513378[28] = 0.0;
   out_3240948688174513378[29] = 0.0;
   out_3240948688174513378[30] = 0.0;
   out_3240948688174513378[31] = 0.0;
   out_3240948688174513378[32] = 0.0;
   out_3240948688174513378[33] = 0.0;
   out_3240948688174513378[34] = 0.0;
   out_3240948688174513378[35] = 0.0;
   out_3240948688174513378[36] = 1.0;
   out_3240948688174513378[37] = 0.0;
   out_3240948688174513378[38] = 0.0;
   out_3240948688174513378[39] = 0.0;
   out_3240948688174513378[40] = 0.0;
   out_3240948688174513378[41] = 0.0;
   out_3240948688174513378[42] = 0.0;
   out_3240948688174513378[43] = 0.0;
   out_3240948688174513378[44] = 0.0;
   out_3240948688174513378[45] = 1.0;
   out_3240948688174513378[46] = 0.0;
   out_3240948688174513378[47] = 0.0;
   out_3240948688174513378[48] = 0.0;
   out_3240948688174513378[49] = 0.0;
   out_3240948688174513378[50] = 0.0;
   out_3240948688174513378[51] = 0.0;
   out_3240948688174513378[52] = 0.0;
   out_3240948688174513378[53] = 0.0;
   out_3240948688174513378[54] = 1.0;
   out_3240948688174513378[55] = 0.0;
   out_3240948688174513378[56] = 0.0;
   out_3240948688174513378[57] = 0.0;
   out_3240948688174513378[58] = 0.0;
   out_3240948688174513378[59] = 0.0;
   out_3240948688174513378[60] = 0.0;
   out_3240948688174513378[61] = 0.0;
   out_3240948688174513378[62] = 0.0;
   out_3240948688174513378[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_2461665730336482261) {
   out_2461665730336482261[0] = state[0];
   out_2461665730336482261[1] = state[1];
   out_2461665730336482261[2] = state[2];
   out_2461665730336482261[3] = state[3];
   out_2461665730336482261[4] = state[4];
   out_2461665730336482261[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2461665730336482261[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2461665730336482261[7] = state[7];
}
void F_fun(double *state, double dt, double *out_7236251358662455363) {
   out_7236251358662455363[0] = 1;
   out_7236251358662455363[1] = 0;
   out_7236251358662455363[2] = 0;
   out_7236251358662455363[3] = 0;
   out_7236251358662455363[4] = 0;
   out_7236251358662455363[5] = 0;
   out_7236251358662455363[6] = 0;
   out_7236251358662455363[7] = 0;
   out_7236251358662455363[8] = 0;
   out_7236251358662455363[9] = 1;
   out_7236251358662455363[10] = 0;
   out_7236251358662455363[11] = 0;
   out_7236251358662455363[12] = 0;
   out_7236251358662455363[13] = 0;
   out_7236251358662455363[14] = 0;
   out_7236251358662455363[15] = 0;
   out_7236251358662455363[16] = 0;
   out_7236251358662455363[17] = 0;
   out_7236251358662455363[18] = 1;
   out_7236251358662455363[19] = 0;
   out_7236251358662455363[20] = 0;
   out_7236251358662455363[21] = 0;
   out_7236251358662455363[22] = 0;
   out_7236251358662455363[23] = 0;
   out_7236251358662455363[24] = 0;
   out_7236251358662455363[25] = 0;
   out_7236251358662455363[26] = 0;
   out_7236251358662455363[27] = 1;
   out_7236251358662455363[28] = 0;
   out_7236251358662455363[29] = 0;
   out_7236251358662455363[30] = 0;
   out_7236251358662455363[31] = 0;
   out_7236251358662455363[32] = 0;
   out_7236251358662455363[33] = 0;
   out_7236251358662455363[34] = 0;
   out_7236251358662455363[35] = 0;
   out_7236251358662455363[36] = 1;
   out_7236251358662455363[37] = 0;
   out_7236251358662455363[38] = 0;
   out_7236251358662455363[39] = 0;
   out_7236251358662455363[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7236251358662455363[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7236251358662455363[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7236251358662455363[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7236251358662455363[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7236251358662455363[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7236251358662455363[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7236251358662455363[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7236251358662455363[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7236251358662455363[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7236251358662455363[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7236251358662455363[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7236251358662455363[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7236251358662455363[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7236251358662455363[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7236251358662455363[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7236251358662455363[56] = 0;
   out_7236251358662455363[57] = 0;
   out_7236251358662455363[58] = 0;
   out_7236251358662455363[59] = 0;
   out_7236251358662455363[60] = 0;
   out_7236251358662455363[61] = 0;
   out_7236251358662455363[62] = 0;
   out_7236251358662455363[63] = 1;
}
void h_25(double *state, double *unused, double *out_7901409450167579230) {
   out_7901409450167579230[0] = state[6];
}
void H_25(double *state, double *unused, double *out_8923886185669432917) {
   out_8923886185669432917[0] = 0;
   out_8923886185669432917[1] = 0;
   out_8923886185669432917[2] = 0;
   out_8923886185669432917[3] = 0;
   out_8923886185669432917[4] = 0;
   out_8923886185669432917[5] = 0;
   out_8923886185669432917[6] = 1;
   out_8923886185669432917[7] = 0;
}
void h_24(double *state, double *unused, double *out_6400597394679481289) {
   out_6400597394679481289[0] = state[4];
   out_6400597394679481289[1] = state[5];
}
void H_24(double *state, double *unused, double *out_7661497121473692785) {
   out_7661497121473692785[0] = 0;
   out_7661497121473692785[1] = 0;
   out_7661497121473692785[2] = 0;
   out_7661497121473692785[3] = 0;
   out_7661497121473692785[4] = 1;
   out_7661497121473692785[5] = 0;
   out_7661497121473692785[6] = 0;
   out_7661497121473692785[7] = 0;
   out_7661497121473692785[8] = 0;
   out_7661497121473692785[9] = 0;
   out_7661497121473692785[10] = 0;
   out_7661497121473692785[11] = 0;
   out_7661497121473692785[12] = 0;
   out_7661497121473692785[13] = 1;
   out_7661497121473692785[14] = 0;
   out_7661497121473692785[15] = 0;
}
void h_30(double *state, double *unused, double *out_3084775859293084061) {
   out_3084775859293084061[0] = state[4];
}
void H_30(double *state, double *unused, double *out_91041022909064581) {
   out_91041022909064581[0] = 0;
   out_91041022909064581[1] = 0;
   out_91041022909064581[2] = 0;
   out_91041022909064581[3] = 0;
   out_91041022909064581[4] = 1;
   out_91041022909064581[5] = 0;
   out_91041022909064581[6] = 0;
   out_91041022909064581[7] = 0;
}
void h_26(double *state, double *unused, double *out_4028285712194746496) {
   out_4028285712194746496[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3219293030730223644) {
   out_3219293030730223644[0] = 0;
   out_3219293030730223644[1] = 0;
   out_3219293030730223644[2] = 0;
   out_3219293030730223644[3] = 0;
   out_3219293030730223644[4] = 0;
   out_3219293030730223644[5] = 0;
   out_3219293030730223644[6] = 0;
   out_3219293030730223644[7] = 1;
}
void h_27(double *state, double *unused, double *out_2899736182506733585) {
   out_2899736182506733585[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1378623010745689893) {
   out_1378623010745689893[0] = 0;
   out_1378623010745689893[1] = 0;
   out_1378623010745689893[2] = 0;
   out_1378623010745689893[3] = 1;
   out_1378623010745689893[4] = 0;
   out_1378623010745689893[5] = 0;
   out_1378623010745689893[6] = 0;
   out_1378623010745689893[7] = 0;
}
void h_29(double *state, double *unused, double *out_4350367575135616752) {
   out_4350367575135616752[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3131189496596960291) {
   out_3131189496596960291[0] = 0;
   out_3131189496596960291[1] = 1;
   out_3131189496596960291[2] = 0;
   out_3131189496596960291[3] = 0;
   out_3131189496596960291[4] = 0;
   out_3131189496596960291[5] = 0;
   out_3131189496596960291[6] = 0;
   out_3131189496596960291[7] = 0;
}
void h_28(double *state, double *unused, double *out_6660478629896449209) {
   out_6660478629896449209[0] = state[5];
   out_6660478629896449209[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6402364645397658906) {
   out_6402364645397658906[0] = 0;
   out_6402364645397658906[1] = 0;
   out_6402364645397658906[2] = 0;
   out_6402364645397658906[3] = 0;
   out_6402364645397658906[4] = 0;
   out_6402364645397658906[5] = 1;
   out_6402364645397658906[6] = 0;
   out_6402364645397658906[7] = 0;
   out_6402364645397658906[8] = 0;
   out_6402364645397658906[9] = 0;
   out_6402364645397658906[10] = 0;
   out_6402364645397658906[11] = 0;
   out_6402364645397658906[12] = 0;
   out_6402364645397658906[13] = 0;
   out_6402364645397658906[14] = 1;
   out_6402364645397658906[15] = 0;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_3837833141682255610) {
  err_fun(nom_x, delta_x, out_3837833141682255610);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3068960279000246371) {
  inv_err_fun(nom_x, true_x, out_3068960279000246371);
}
void car_H_mod_fun(double *state, double *out_3240948688174513378) {
  H_mod_fun(state, out_3240948688174513378);
}
void car_f_fun(double *state, double dt, double *out_2461665730336482261) {
  f_fun(state,  dt, out_2461665730336482261);
}
void car_F_fun(double *state, double dt, double *out_7236251358662455363) {
  F_fun(state,  dt, out_7236251358662455363);
}
void car_h_25(double *state, double *unused, double *out_7901409450167579230) {
  h_25(state, unused, out_7901409450167579230);
}
void car_H_25(double *state, double *unused, double *out_8923886185669432917) {
  H_25(state, unused, out_8923886185669432917);
}
void car_h_24(double *state, double *unused, double *out_6400597394679481289) {
  h_24(state, unused, out_6400597394679481289);
}
void car_H_24(double *state, double *unused, double *out_7661497121473692785) {
  H_24(state, unused, out_7661497121473692785);
}
void car_h_30(double *state, double *unused, double *out_3084775859293084061) {
  h_30(state, unused, out_3084775859293084061);
}
void car_H_30(double *state, double *unused, double *out_91041022909064581) {
  H_30(state, unused, out_91041022909064581);
}
void car_h_26(double *state, double *unused, double *out_4028285712194746496) {
  h_26(state, unused, out_4028285712194746496);
}
void car_H_26(double *state, double *unused, double *out_3219293030730223644) {
  H_26(state, unused, out_3219293030730223644);
}
void car_h_27(double *state, double *unused, double *out_2899736182506733585) {
  h_27(state, unused, out_2899736182506733585);
}
void car_H_27(double *state, double *unused, double *out_1378623010745689893) {
  H_27(state, unused, out_1378623010745689893);
}
void car_h_29(double *state, double *unused, double *out_4350367575135616752) {
  h_29(state, unused, out_4350367575135616752);
}
void car_H_29(double *state, double *unused, double *out_3131189496596960291) {
  H_29(state, unused, out_3131189496596960291);
}
void car_h_28(double *state, double *unused, double *out_6660478629896449209) {
  h_28(state, unused, out_6660478629896449209);
}
void car_H_28(double *state, double *unused, double *out_6402364645397658906) {
  H_28(state, unused, out_6402364645397658906);
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
