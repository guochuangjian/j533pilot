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
void car_err_fun(double *nom_x, double *delta_x, double *out_5826554309974982227);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_1829394537637158322);
void car_H_mod_fun(double *state, double *out_7631383496007094190);
void car_f_fun(double *state, double dt, double *out_6042263353151296703);
void car_F_fun(double *state, double dt, double *out_770756480427470884);
void car_h_25(double *state, double *unused, double *out_5361102727716958993);
void car_H_25(double *state, double *unused, double *out_3207636012344091342);
void car_h_24(double *state, double *unused, double *out_2493276794756230534);
void car_H_24(double *state, double *unused, double *out_3137807261716376900);
void car_h_30(double *state, double *unused, double *out_4902901365558397331);
void car_H_30(double *state, double *unused, double *out_6406262898605091938);
void car_h_26(double *state, double *unused, double *out_3893408098707599809);
void car_H_26(double *state, double *unused, double *out_8912229167283300615);
void car_h_27(double *state, double *unused, double *out_5964418612393436510);
void car_H_27(double *state, double *unused, double *out_7693844886441717250);
void car_h_29(double *state, double *unused, double *out_1866115314114064646);
void car_H_29(double *state, double *unused, double *out_5831704284749555763);
void car_h_28(double *state, double *unused, double *out_3369080615822404081);
void car_H_28(double *state, double *unused, double *out_8291104333756679850);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}