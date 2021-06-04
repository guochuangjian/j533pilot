#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_3(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_19(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_1407597752604516151);
void live_err_fun(double *nom_x, double *delta_x, double *out_2755351900306618068);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8740450469725871824);
void live_H_mod_fun(double *state, double *out_3507290464964589302);
void live_f_fun(double *state, double dt, double *out_4045014662031557180);
void live_F_fun(double *state, double dt, double *out_7030749745690010680);
void live_h_3(double *state, double *unused, double *out_1579097727061021667);
void live_H_3(double *state, double *unused, double *out_7023550732714914920);
void live_h_4(double *state, double *unused, double *out_9097939223678211760);
void live_H_4(double *state, double *unused, double *out_418623243602265916);
void live_h_9(double *state, double *unused, double *out_7191140075104019335);
void live_H_9(double *state, double *unused, double *out_8659585997601572905);
void live_h_10(double *state, double *unused, double *out_58680644405209996);
void live_H_10(double *state, double *unused, double *out_7212978363181974909);
void live_h_12(double *state, double *unused, double *out_9143016397561043394);
void live_H_12(double *state, double *unused, double *out_2525178755341475485);
void live_h_31(double *state, double *unused, double *out_2446152173738630444);
void live_H_31(double *state, double *unused, double *out_149911019801944336);
void live_h_32(double *state, double *unused, double *out_8758063247912175626);
void live_H_32(double *state, double *unused, double *out_8259950316454613854);
void live_h_13(double *state, double *unused, double *out_3880516526716149206);
void live_H_13(double *state, double *unused, double *out_7066692656685020581);
void live_h_14(double *state, double *unused, double *out_7191140075104019335);
void live_H_14(double *state, double *unused, double *out_8659585997601572905);
void live_h_19(double *state, double *unused, double *out_4456036503049278004);
void live_H_19(double *state, double *unused, double *out_4226906146499141724);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}