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
void err_fun(double *nom_x, double *delta_x, double *out_2512610070777309265) {
   out_2512610070777309265[0] = delta_x[0] + nom_x[0];
   out_2512610070777309265[1] = delta_x[1] + nom_x[1];
   out_2512610070777309265[2] = delta_x[2] + nom_x[2];
   out_2512610070777309265[3] = delta_x[3] + nom_x[3];
   out_2512610070777309265[4] = delta_x[4] + nom_x[4];
   out_2512610070777309265[5] = delta_x[5] + nom_x[5];
   out_2512610070777309265[6] = delta_x[6] + nom_x[6];
   out_2512610070777309265[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5488762515128777546) {
   out_5488762515128777546[0] = -nom_x[0] + true_x[0];
   out_5488762515128777546[1] = -nom_x[1] + true_x[1];
   out_5488762515128777546[2] = -nom_x[2] + true_x[2];
   out_5488762515128777546[3] = -nom_x[3] + true_x[3];
   out_5488762515128777546[4] = -nom_x[4] + true_x[4];
   out_5488762515128777546[5] = -nom_x[5] + true_x[5];
   out_5488762515128777546[6] = -nom_x[6] + true_x[6];
   out_5488762515128777546[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_6473749384078044496) {
   out_6473749384078044496[0] = 1.0;
   out_6473749384078044496[1] = 0.0;
   out_6473749384078044496[2] = 0.0;
   out_6473749384078044496[3] = 0.0;
   out_6473749384078044496[4] = 0.0;
   out_6473749384078044496[5] = 0.0;
   out_6473749384078044496[6] = 0.0;
   out_6473749384078044496[7] = 0.0;
   out_6473749384078044496[8] = 0.0;
   out_6473749384078044496[9] = 1.0;
   out_6473749384078044496[10] = 0.0;
   out_6473749384078044496[11] = 0.0;
   out_6473749384078044496[12] = 0.0;
   out_6473749384078044496[13] = 0.0;
   out_6473749384078044496[14] = 0.0;
   out_6473749384078044496[15] = 0.0;
   out_6473749384078044496[16] = 0.0;
   out_6473749384078044496[17] = 0.0;
   out_6473749384078044496[18] = 1.0;
   out_6473749384078044496[19] = 0.0;
   out_6473749384078044496[20] = 0.0;
   out_6473749384078044496[21] = 0.0;
   out_6473749384078044496[22] = 0.0;
   out_6473749384078044496[23] = 0.0;
   out_6473749384078044496[24] = 0.0;
   out_6473749384078044496[25] = 0.0;
   out_6473749384078044496[26] = 0.0;
   out_6473749384078044496[27] = 1.0;
   out_6473749384078044496[28] = 0.0;
   out_6473749384078044496[29] = 0.0;
   out_6473749384078044496[30] = 0.0;
   out_6473749384078044496[31] = 0.0;
   out_6473749384078044496[32] = 0.0;
   out_6473749384078044496[33] = 0.0;
   out_6473749384078044496[34] = 0.0;
   out_6473749384078044496[35] = 0.0;
   out_6473749384078044496[36] = 1.0;
   out_6473749384078044496[37] = 0.0;
   out_6473749384078044496[38] = 0.0;
   out_6473749384078044496[39] = 0.0;
   out_6473749384078044496[40] = 0.0;
   out_6473749384078044496[41] = 0.0;
   out_6473749384078044496[42] = 0.0;
   out_6473749384078044496[43] = 0.0;
   out_6473749384078044496[44] = 0.0;
   out_6473749384078044496[45] = 1.0;
   out_6473749384078044496[46] = 0.0;
   out_6473749384078044496[47] = 0.0;
   out_6473749384078044496[48] = 0.0;
   out_6473749384078044496[49] = 0.0;
   out_6473749384078044496[50] = 0.0;
   out_6473749384078044496[51] = 0.0;
   out_6473749384078044496[52] = 0.0;
   out_6473749384078044496[53] = 0.0;
   out_6473749384078044496[54] = 1.0;
   out_6473749384078044496[55] = 0.0;
   out_6473749384078044496[56] = 0.0;
   out_6473749384078044496[57] = 0.0;
   out_6473749384078044496[58] = 0.0;
   out_6473749384078044496[59] = 0.0;
   out_6473749384078044496[60] = 0.0;
   out_6473749384078044496[61] = 0.0;
   out_6473749384078044496[62] = 0.0;
   out_6473749384078044496[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_2454989366969667135) {
   out_2454989366969667135[0] = state[0];
   out_2454989366969667135[1] = state[1];
   out_2454989366969667135[2] = state[2];
   out_2454989366969667135[3] = state[3];
   out_2454989366969667135[4] = state[4];
   out_2454989366969667135[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2454989366969667135[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2454989366969667135[7] = state[7];
}
void F_fun(double *state, double dt, double *out_6995264623040855883) {
   out_6995264623040855883[0] = 1;
   out_6995264623040855883[1] = 0;
   out_6995264623040855883[2] = 0;
   out_6995264623040855883[3] = 0;
   out_6995264623040855883[4] = 0;
   out_6995264623040855883[5] = 0;
   out_6995264623040855883[6] = 0;
   out_6995264623040855883[7] = 0;
   out_6995264623040855883[8] = 0;
   out_6995264623040855883[9] = 1;
   out_6995264623040855883[10] = 0;
   out_6995264623040855883[11] = 0;
   out_6995264623040855883[12] = 0;
   out_6995264623040855883[13] = 0;
   out_6995264623040855883[14] = 0;
   out_6995264623040855883[15] = 0;
   out_6995264623040855883[16] = 0;
   out_6995264623040855883[17] = 0;
   out_6995264623040855883[18] = 1;
   out_6995264623040855883[19] = 0;
   out_6995264623040855883[20] = 0;
   out_6995264623040855883[21] = 0;
   out_6995264623040855883[22] = 0;
   out_6995264623040855883[23] = 0;
   out_6995264623040855883[24] = 0;
   out_6995264623040855883[25] = 0;
   out_6995264623040855883[26] = 0;
   out_6995264623040855883[27] = 1;
   out_6995264623040855883[28] = 0;
   out_6995264623040855883[29] = 0;
   out_6995264623040855883[30] = 0;
   out_6995264623040855883[31] = 0;
   out_6995264623040855883[32] = 0;
   out_6995264623040855883[33] = 0;
   out_6995264623040855883[34] = 0;
   out_6995264623040855883[35] = 0;
   out_6995264623040855883[36] = 1;
   out_6995264623040855883[37] = 0;
   out_6995264623040855883[38] = 0;
   out_6995264623040855883[39] = 0;
   out_6995264623040855883[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_6995264623040855883[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_6995264623040855883[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6995264623040855883[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_6995264623040855883[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_6995264623040855883[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_6995264623040855883[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_6995264623040855883[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_6995264623040855883[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_6995264623040855883[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_6995264623040855883[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6995264623040855883[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6995264623040855883[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_6995264623040855883[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_6995264623040855883[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_6995264623040855883[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_6995264623040855883[56] = 0;
   out_6995264623040855883[57] = 0;
   out_6995264623040855883[58] = 0;
   out_6995264623040855883[59] = 0;
   out_6995264623040855883[60] = 0;
   out_6995264623040855883[61] = 0;
   out_6995264623040855883[62] = 0;
   out_6995264623040855883[63] = 1;
}
void h_25(double *state, double *unused, double *out_3168366076751619319) {
   out_3168366076751619319[0] = state[6];
}
void H_25(double *state, double *unused, double *out_5848059173098688088) {
   out_5848059173098688088[0] = 0;
   out_5848059173098688088[1] = 0;
   out_5848059173098688088[2] = 0;
   out_5848059173098688088[3] = 0;
   out_5848059173098688088[4] = 0;
   out_5848059173098688088[5] = 0;
   out_5848059173098688088[6] = 1;
   out_5848059173098688088[7] = 0;
}
void h_24(double *state, double *unused, double *out_967367903101647447) {
   out_967367903101647447[0] = state[4];
   out_967367903101647447[1] = state[5];
}
void H_24(double *state, double *unused, double *out_296241696831249030) {
   out_296241696831249030[0] = 0;
   out_296241696831249030[1] = 0;
   out_296241696831249030[2] = 0;
   out_296241696831249030[3] = 0;
   out_296241696831249030[4] = 1;
   out_296241696831249030[5] = 0;
   out_296241696831249030[6] = 0;
   out_296241696831249030[7] = 0;
   out_296241696831249030[8] = 0;
   out_296241696831249030[9] = 0;
   out_296241696831249030[10] = 0;
   out_296241696831249030[11] = 0;
   out_296241696831249030[12] = 0;
   out_296241696831249030[13] = 1;
   out_296241696831249030[14] = 0;
   out_296241696831249030[15] = 0;
}
void h_30(double *state, double *unused, double *out_3022846252530964348) {
   out_3022846252530964348[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1413571393322687880) {
   out_1413571393322687880[0] = 0;
   out_1413571393322687880[1] = 0;
   out_1413571393322687880[2] = 0;
   out_1413571393322687880[3] = 0;
   out_1413571393322687880[4] = 1;
   out_1413571393322687880[5] = 0;
   out_1413571393322687880[6] = 0;
   out_1413571393322687880[7] = 0;
}
void h_26(double *state, double *unused, double *out_8024519520191809993) {
   out_8024519520191809993[0] = state[7];
}
void H_26(double *state, double *unused, double *out_6858891383930847848) {
   out_6858891383930847848[0] = 0;
   out_6858891383930847848[1] = 0;
   out_6858891383930847848[2] = 0;
   out_6858891383930847848[3] = 0;
   out_6858891383930847848[4] = 0;
   out_6858891383930847848[5] = 0;
   out_6858891383930847848[6] = 0;
   out_6858891383930847848[7] = 1;
}
void h_27(double *state, double *unused, double *out_9151289916057761072) {
   out_9151289916057761072[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1697204001825054936) {
   out_1697204001825054936[0] = 0;
   out_1697204001825054936[1] = 0;
   out_1697204001825054936[2] = 0;
   out_1697204001825054936[3] = 1;
   out_1697204001825054936[4] = 0;
   out_1697204001825054936[5] = 0;
   out_1697204001825054936[6] = 0;
   out_1697204001825054936[7] = 0;
}
void h_29(double *state, double *unused, double *out_3769396721100767431) {
   out_3769396721100767431[0] = state[1];
}
void H_29(double *state, double *unused, double *out_839012779467151705) {
   out_839012779467151705[0] = 0;
   out_839012779467151705[1] = 1;
   out_839012779467151705[2] = 0;
   out_839012779467151705[3] = 0;
   out_839012779467151705[4] = 0;
   out_839012779467151705[5] = 0;
   out_839012779467151705[6] = 0;
   out_839012779467151705[7] = 0;
}
void h_28(double *state, double *unused, double *out_5593901768300893881) {
   out_5593901768300893881[0] = state[5];
   out_5593901768300893881[1] = state[6];
}
void H_28(double *state, double *unused, double *out_6721590781405245836) {
   out_6721590781405245836[0] = 0;
   out_6721590781405245836[1] = 0;
   out_6721590781405245836[2] = 0;
   out_6721590781405245836[3] = 0;
   out_6721590781405245836[4] = 0;
   out_6721590781405245836[5] = 1;
   out_6721590781405245836[6] = 0;
   out_6721590781405245836[7] = 0;
   out_6721590781405245836[8] = 0;
   out_6721590781405245836[9] = 0;
   out_6721590781405245836[10] = 0;
   out_6721590781405245836[11] = 0;
   out_6721590781405245836[12] = 0;
   out_6721590781405245836[13] = 0;
   out_6721590781405245836[14] = 1;
   out_6721590781405245836[15] = 0;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_2512610070777309265) {
  err_fun(nom_x, delta_x, out_2512610070777309265);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5488762515128777546) {
  inv_err_fun(nom_x, true_x, out_5488762515128777546);
}
void car_H_mod_fun(double *state, double *out_6473749384078044496) {
  H_mod_fun(state, out_6473749384078044496);
}
void car_f_fun(double *state, double dt, double *out_2454989366969667135) {
  f_fun(state,  dt, out_2454989366969667135);
}
void car_F_fun(double *state, double dt, double *out_6995264623040855883) {
  F_fun(state,  dt, out_6995264623040855883);
}
void car_h_25(double *state, double *unused, double *out_3168366076751619319) {
  h_25(state, unused, out_3168366076751619319);
}
void car_H_25(double *state, double *unused, double *out_5848059173098688088) {
  H_25(state, unused, out_5848059173098688088);
}
void car_h_24(double *state, double *unused, double *out_967367903101647447) {
  h_24(state, unused, out_967367903101647447);
}
void car_H_24(double *state, double *unused, double *out_296241696831249030) {
  H_24(state, unused, out_296241696831249030);
}
void car_h_30(double *state, double *unused, double *out_3022846252530964348) {
  h_30(state, unused, out_3022846252530964348);
}
void car_H_30(double *state, double *unused, double *out_1413571393322687880) {
  H_30(state, unused, out_1413571393322687880);
}
void car_h_26(double *state, double *unused, double *out_8024519520191809993) {
  h_26(state, unused, out_8024519520191809993);
}
void car_H_26(double *state, double *unused, double *out_6858891383930847848) {
  H_26(state, unused, out_6858891383930847848);
}
void car_h_27(double *state, double *unused, double *out_9151289916057761072) {
  h_27(state, unused, out_9151289916057761072);
}
void car_H_27(double *state, double *unused, double *out_1697204001825054936) {
  H_27(state, unused, out_1697204001825054936);
}
void car_h_29(double *state, double *unused, double *out_3769396721100767431) {
  h_29(state, unused, out_3769396721100767431);
}
void car_H_29(double *state, double *unused, double *out_839012779467151705) {
  H_29(state, unused, out_839012779467151705);
}
void car_h_28(double *state, double *unused, double *out_5593901768300893881) {
  h_28(state, unused, out_5593901768300893881);
}
void car_H_28(double *state, double *unused, double *out_6721590781405245836) {
  H_28(state, unused, out_6721590781405245836);
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
