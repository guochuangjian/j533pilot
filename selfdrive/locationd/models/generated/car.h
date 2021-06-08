#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3837833141682255610);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_3068960279000246371);
void car_H_mod_fun(double *state, double *out_3240948688174513378);
void car_f_fun(double *state, double dt, double *out_2461665730336482261);
void car_F_fun(double *state, double dt, double *out_7236251358662455363);
void car_h_25(double *state, double *unused, double *out_7901409450167579230);
void car_H_25(double *state, double *unused, double *out_8923886185669432917);
void car_h_24(double *state, double *unused, double *out_6400597394679481289);
void car_H_24(double *state, double *unused, double *out_7661497121473692785);
void car_h_30(double *state, double *unused, double *out_3084775859293084061);
void car_H_30(double *state, double *unused, double *out_91041022909064581);
void car_h_26(double *state, double *unused, double *out_4028285712194746496);
void car_H_26(double *state, double *unused, double *out_3219293030730223644);
void car_h_27(double *state, double *unused, double *out_2899736182506733585);
void car_H_27(double *state, double *unused, double *out_1378623010745689893);
void car_h_29(double *state, double *unused, double *out_4350367575135616752);
void car_H_29(double *state, double *unused, double *out_3131189496596960291);
void car_h_28(double *state, double *unused, double *out_6660478629896449209);
void car_H_28(double *state, double *unused, double *out_6402364645397658906);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}