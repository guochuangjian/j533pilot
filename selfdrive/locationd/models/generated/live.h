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
void live_H(double *in_vec, double *out_3112248656344267024);
void live_err_fun(double *nom_x, double *delta_x, double *out_6514478479377154554);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6562831803697931829);
void live_H_mod_fun(double *state, double *out_3485285672998159801);
void live_f_fun(double *state, double dt, double *out_2252742358831627050);
void live_F_fun(double *state, double dt, double *out_1390635145567685639);
void live_h_3(double *state, double *unused, double *out_8778689902229900358);
void live_H_3(double *state, double *unused, double *out_7071063248221903837);
void live_h_4(double *state, double *unused, double *out_7412584150661195934);
void live_H_4(double *state, double *unused, double *out_8513600165493837225);
void live_h_9(double *state, double *unused, double *out_2165147619284547662);
void live_H_9(double *state, double *unused, double *out_854934667011875570);
void live_h_10(double *state, double *unused, double *out_6437534819363874302);
void live_H_10(double *state, double *unused, double *out_4283058499254228570);
void live_h_12(double *state, double *unused, double *out_294143016965715036);
void live_H_12(double *state, double *unused, double *out_6407044653754627656);
void live_h_31(double *state, double *unused, double *out_2742120223357459888);
void live_H_31(double *state, double *unused, double *out_9082134428898047477);
void live_h_32(double *state, double *unused, double *out_2553674345415990180);
void live_H_32(double *state, double *unused, double *out_1683532725355929790);
void live_h_13(double *state, double *unused, double *out_2051282062101304248);
void live_H_13(double *state, double *unused, double *out_6082984043611976909);
void live_h_14(double *state, double *unused, double *out_2165147619284547662);
void live_H_14(double *state, double *unused, double *out_854934667011875570);
void live_h_19(double *state, double *unused, double *out_2529107655439357541);
void live_H_19(double *state, double *unused, double *out_5287614518114306751);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}