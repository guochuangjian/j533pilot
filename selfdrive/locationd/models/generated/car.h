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
void car_err_fun(double *nom_x, double *delta_x, double *out_2512610070777309265);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5488762515128777546);
void car_H_mod_fun(double *state, double *out_6473749384078044496);
void car_f_fun(double *state, double dt, double *out_2454989366969667135);
void car_F_fun(double *state, double dt, double *out_6995264623040855883);
void car_h_25(double *state, double *unused, double *out_3168366076751619319);
void car_H_25(double *state, double *unused, double *out_5848059173098688088);
void car_h_24(double *state, double *unused, double *out_967367903101647447);
void car_H_24(double *state, double *unused, double *out_296241696831249030);
void car_h_30(double *state, double *unused, double *out_3022846252530964348);
void car_H_30(double *state, double *unused, double *out_1413571393322687880);
void car_h_26(double *state, double *unused, double *out_8024519520191809993);
void car_H_26(double *state, double *unused, double *out_6858891383930847848);
void car_h_27(double *state, double *unused, double *out_9151289916057761072);
void car_H_27(double *state, double *unused, double *out_1697204001825054936);
void car_h_29(double *state, double *unused, double *out_3769396721100767431);
void car_H_29(double *state, double *unused, double *out_839012779467151705);
void car_h_28(double *state, double *unused, double *out_5593901768300893881);
void car_H_28(double *state, double *unused, double *out_6721590781405245836);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}