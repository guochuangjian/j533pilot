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
void live_H(double *in_vec, double *out_1972773286786522443);
void live_err_fun(double *nom_x, double *delta_x, double *out_459408300567769233);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_246932384933284732);
void live_H_mod_fun(double *state, double *out_2832780683553674837);
void live_f_fun(double *state, double dt, double *out_8655927765837852660);
void live_F_fun(double *state, double dt, double *out_3108770954200234939);
void live_h_3(double *state, double *unused, double *out_4082749196527027513);
void live_H_3(double *state, double *unused, double *out_1912653736884298123);
void live_h_4(double *state, double *unused, double *out_2290324257456542897);
void live_H_4(double *state, double *unused, double *out_1656412342597990744);
void live_h_9(double *state, double *unused, double *out_1515964259684499922);
void live_H_9(double *state, double *unused, double *out_7421796898605848077);
void live_h_10(double *state, double *unused, double *out_7542071773270820038);
void live_H_10(double *state, double *unused, double *out_7505733408038929869);
void live_h_12(double *state, double *unused, double *out_3438056662012513987);
void live_H_12(double *state, double *unused, double *out_3762967854337200313);
void live_h_31(double *state, double *unused, double *out_5155955351434917179);
void live_H_31(double *state, double *unused, double *out_1087878079193780492);
void live_h_32(double *state, double *unused, double *out_2495499613004245457);
void live_H_32(double *state, double *unused, double *out_5212487665066620694);
void live_h_13(double *state, double *unused, double *out_8487153226629105807);
void live_H_13(double *state, double *unused, double *out_2623981466124861801);
void live_h_14(double *state, double *unused, double *out_1515964259684499922);
void live_H_14(double *state, double *unused, double *out_7421796898605848077);
void live_h_19(double *state, double *unused, double *out_980967387768816414);
void live_H_19(double *state, double *unused, double *out_2989117047503416896);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}