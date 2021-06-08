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
void live_H(double *in_vec, double *out_856236467480591779);
void live_err_fun(double *nom_x, double *delta_x, double *out_6469533078858542841);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_9185883606417194634);
void live_H_mod_fun(double *state, double *out_4511572670331821600);
void live_f_fun(double *state, double dt, double *out_3241693131083666832);
void live_F_fun(double *state, double dt, double *out_4642879992605947154);
void live_h_3(double *state, double *unused, double *out_4714584368116173582);
void live_H_3(double *state, double *unused, double *out_7920832901572434099);
void live_h_4(double *state, double *unused, double *out_2288116132464842531);
void live_H_4(double *state, double *unused, double *out_4804106790858952924);
void live_h_9(double *state, double *unused, double *out_8439819811893802402);
void live_H_9(double *state, double *unused, double *out_2481601246988096954);
void live_h_10(double *state, double *unused, double *out_4486324826627784159);
void live_H_10(double *state, double *unused, double *out_2939160307098811205);
void live_h_12(double *state, double *unused, double *out_2578323031660001065);
void live_H_12(double *state, double *unused, double *out_2697551279119743355);
void live_h_31(double *state, double *unused, double *out_3372423484097720964);
void live_H_31(double *state, double *unused, double *out_5372641054263163176);
void live_h_32(double *state, double *unused, double *out_4209223614226194236);
void live_H_32(double *state, double *unused, double *out_6773737275185987254);
void live_h_13(double *state, double *unused, double *out_4795803117852923289);
void live_H_13(double *state, double *unused, double *out_3484022974191141879);
void live_h_14(double *state, double *unused, double *out_8439819811893802402);
void live_H_14(double *state, double *unused, double *out_2481601246988096954);
void live_h_19(double *state, double *unused, double *out_3231881463738757134);
void live_H_19(double *state, double *unused, double *out_8997107892749191052);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}