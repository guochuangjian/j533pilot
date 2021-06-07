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
void car_err_fun(double *nom_x, double *delta_x, double *out_1988559408514834698);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5716367333266523332);
void car_H_mod_fun(double *state, double *out_7495240528319181325);
void car_f_fun(double *state, double dt, double *out_3177971419004634205);
void car_F_fun(double *state, double dt, double *out_2211924309575664307);
void car_h_25(double *state, double *unused, double *out_2083586975813753818);
void car_H_25(double *state, double *unused, double *out_3475252786579172581);
void car_h_24(double *state, double *unused, double *out_5688669739645459629);
void car_H_24(double *state, double *unused, double *out_3545081537206887023);
void car_h_30(double *state, double *unused, double *out_1482168818200390264);
void car_H_30(double *state, double *unused, double *out_5357592376181195755);
void car_h_26(double *state, double *unused, double *out_878852933523912747);
void car_H_26(double *state, double *unused, double *out_2229340368360036692);
void car_h_27(double *state, double *unused, double *out_4791496524062440986);
void car_H_27(double *state, double *unused, double *out_4070010388344570443);
void car_h_29(double *state, double *unused, double *out_7499800353806450501);
void car_H_29(double *state, double *unused, double *out_5932150990036731930);
void car_h_28(double *state, double *unused, double *out_1545142520145377979);
void car_H_28(double *state, double *unused, double *out_3472750941029607843);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}