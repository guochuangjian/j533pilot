#include "selfdrive/ui/paint_extend.h"
#include "selfdrive/ui/ui.h"
#include <math.h>

void ui_draw_speed_image(const UIState *s, int x, int y, int w, int h, const char *name) {
  nvgBeginPath(s->vg);
  NVGpaint imgPaint = nvgImagePattern(s->vg, x, y, w, h, 0, s->images.at(name), 1.0f);
  nvgRect(s->vg, x, y, w, h);
  nvgFillPaint(s->vg, imgPaint);
  nvgFill(s->vg);
}
//===== draw text =====
void ui_draw_hud_text(UIState *s,
                        const int x,
                        const int y,
                        const char* font_string,
                        const int font_size,
                        const NVGcolor font_color) {
  nvgFontSize(s->vg, font_size);
  nvgFillColor(s->vg, font_color);
  nvgTextAlign(s->vg, NVG_ALIGN_LEFT|NVG_ALIGN_TOP);
  nvgText(s->vg, x, y, font_string, NULL);
}
//===== draw box/title/value/unit =====
void ui_draw_hud_box(UIState *s,
                    const int x,
                    const int y,
                    const int w,
                    const int h) {
  nvgBeginPath(s->vg);
  nvgRoundedRect(s->vg, x, y, h, w, 20);
  nvgStrokeColor(s->vg, COLOR_WHITE_ALPHA(80));
  nvgStrokeWidth(s->vg, 6);
  nvgStroke(s->vg);
}
void ui_draw_hud_title(UIState *s,
                        const int x,
                        const int y,
                        const int w,
                        const int h,
                        const char* font_string,
                        const int font_size,
                        const NVGcolor font_color,
                        const NVGalign font_algin) {
  nvgFontSize(s->vg, font_size);
  nvgFillColor(s->vg, font_color);
  nvgTextAlign(s->vg, font_algin);
  nvgText(s->vg, x + w / 2, y + 40, font_string, NULL);
}
void ui_draw_hud_value(UIState *s,
                        const int x,
                        const int y,
                        const int w,
                        const int h,
                        const char* font_string,
                        const int font_size,
                        const NVGcolor font_color,
                        const NVGalign font_algin) {
  nvgFontSize(s->vg, font_size);
  nvgFillColor(s->vg, font_color);
  nvgTextAlign(s->vg, font_algin);
  nvgText(s->vg, x + w / 2, y + h - 70, font_string, NULL);
}

void ui_draw_hud_unit(UIState *s,
                        const int x,
                        const int y,
                        const int w,
                        const int h,
                        const char* font_string,
                        const int font_size,
                        const NVGcolor font_color,
                        const NVGalign font_algin) {
  nvgFontSize(s->vg, font_size);
  nvgFillColor(s->vg, font_color);
  nvgTextAlign(s->vg, font_algin);
  nvgText(s->vg, x + w / 2, y + h - 20, font_string, NULL);
}

void ui_draw_hud_infobox(UIState *s,
                        const int x,
                        const int y,
                        const int w,
                        const int h,
                        const char* title_font_string,
                        const int title_font_size,
                        const NVGcolor title_font_color,
                        const NVGalign title_font_algin,
                        const char* value_font_string,
                        const int value_font_size,
                        const NVGcolor value_font_color,
                        const NVGalign value_font_algin,
                        const char* unit_font_string,
                        const int unit_font_size,
                        const NVGcolor unit_font_color,
                        const NVGalign unit_font_algin) {
  ui_draw_hud_box(s, x, y, w, h);
  ui_draw_hud_title(s, x, y, w, h, title_font_string, title_font_size, title_font_color, title_font_algin);
  ui_draw_hud_value(s, x, y, w, h, value_font_string, value_font_size, value_font_color, value_font_algin);
  ui_draw_hud_unit(s, x, y, w, h, unit_font_string, unit_font_size, unit_font_color, unit_font_algin);
}

//===== draw top hud 1/2/3/4 =====
void ui_draw_top_hud_infobox1(UIState *s) {
}

void ui_draw_top_hud_infobox2(UIState *s) {
}

void ui_draw_top_hud_infobox3(UIState *s) {
}

void ui_draw_top_hud_infobox4(UIState *s) {
}


//===== draw hud top/bottom/left/right
void ui_draw_top_hud(UIState *s) {

  ui_draw_top_hud_infobox1(s);
  ui_draw_top_hud_infobox2(s);
  ui_draw_top_hud_infobox3(s);
  ui_draw_top_hud_infobox4(s);
}

void ui_draw_bottom_hud(UIState *s) {
}

void ui_draw_left_hud_infobox1(UIState *s) {
  int sidebar_fit_x = 0;
  char value[16];
  NVGcolor value_font_color = COLOR_WHITE_ALPHA(200);
  float steeringAngleDeg = s->scene.car_state.getSteeringAngleDeg();

  //Fit sidebar screen
  sidebar_fit_x = s->viz_rect.x + hud_left_1_x;

  if(((int)(steeringAngleDeg) < -30) || ((int)(steeringAngleDeg) > 30)) {
    value_font_color = nvgRGBA(0, 255, 255, 200);
  }
  if(((int)(steeringAngleDeg) < -60) || ((int)(steeringAngleDeg) > 60)) {
    value_font_color = nvgRGBA(0, 255, 0, 200);
  }
  if(((int)(steeringAngleDeg) < -90) || ((int)(steeringAngleDeg) > 90)) {
    value_font_color = nvgRGBA(255, 255, 0, 200);
  }
  if(((int)(steeringAngleDeg) < -120) || ((int)(steeringAngleDeg) > 120)) {
    value_font_color = nvgRGBA(255, 127, 0, 200);
  }
  if(((int)(steeringAngleDeg) < -150) || ((int)(steeringAngleDeg) > 150)) {
    value_font_color = nvgRGBA(255, 0, 255, 200);
  }
  if(((int)(steeringAngleDeg) < -180) || ((int)(steeringAngleDeg) > 180)) {
    value_font_color = nvgRGBA(255, 0, 0, 200);
  }

  //value
  snprintf(value, sizeof(value), "%.2f°",(steeringAngleDeg));
  if(s->scene.car_state.getSteeringPressed()) {
    ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_1_y, hud_left_1_w, hud_left_1_h,
                          "User Angle", 40, COLOR_YELLOW, NVG_ALIGN_CENTER,
                          value, 60, value_font_color, NVG_ALIGN_CENTER,
                          "degree", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
  } else if(s->scene.car_control.getAvailableFulltimeLka()) {
    ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_1_y, hud_left_1_w, hud_left_1_h,
                          "OP Angle", 48, COLOR_GREEN, NVG_ALIGN_CENTER,
                          value, 60, value_font_color, NVG_ALIGN_CENTER,
                          "degree", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
  } else {
    ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_1_y, hud_left_1_w, hud_left_1_h,
                          "Angle", 48, COLOR_WHITE, NVG_ALIGN_CENTER,
                          value, 60, value_font_color, NVG_ALIGN_CENTER,
                          "degree", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
  }
}

void ui_draw_left_hud_infobox2(UIState *s) {
  int sidebar_fit_x = 0;
  char value[16];

  //Fit sidebar screen
  sidebar_fit_x = s->viz_rect.x + hud_left_2_x;

  if(s->scene.car_state.getBrake()>0) {
      snprintf(value, sizeof(value), "%.2f°", s->scene.car_state.getBrake());
      ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_2_y, hud_left_2_w, hud_left_2_h,
                        "User Brake", 40, COLOR_YELLOW, NVG_ALIGN_CENTER,
                        value, 72, COLOR_RED, NVG_ALIGN_CENTER,
                        "", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
  } else {
    if(s->scene.car_state.getGas()>0) {
      snprintf(value, sizeof(value), "%.2f°",(s->scene.car_state.getGas()));
      if(s->scene.car_state.getGasPressed()) {
          ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_2_y, hud_left_2_w, hud_left_2_h,
                            "User Gas", 40, COLOR_YELLOW, NVG_ALIGN_CENTER,
                            value, 72, COLOR_GREEN, NVG_ALIGN_CENTER,
                            "", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
      } else {
          ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_2_y, hud_left_2_w, hud_left_2_h,
                            "OP Gas", 40, COLOR_GREEN, NVG_ALIGN_CENTER,
                            value, 72, COLOR_GREEN, NVG_ALIGN_CENTER,
                            "", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
      }
    } else {
      ui_draw_hud_infobox(s, sidebar_fit_x, hud_left_2_y, hud_left_2_w, hud_left_2_h,
                            "Gas/Brake", 40, COLOR_WHITE, NVG_ALIGN_CENTER,
                            "", 72, COLOR_WHITE, NVG_ALIGN_CENTER,
                            "", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
    }
  }
}

void ui_draw_left_hud(UIState *s) {
  ui_draw_left_hud_infobox1(s);
  ui_draw_left_hud_infobox2(s);
}

void ui_draw_right_hud_infobox1(UIState *s) {
  float d_rel = s->scene.lead_data[0].getDRel();
  //float v_rel = s->scene.lead_data[0].getVRel();
  char value[16];

  snprintf(value, sizeof(value), "%.2f", d_rel);
  ui_draw_hud_infobox(s, hud_right_1_x, hud_right_1_y, hud_right_1_w, hud_right_1_h,
                            "Lead Dis", 40, COLOR_WHITE, NVG_ALIGN_CENTER,
                            "", 72, COLOR_WHITE, NVG_ALIGN_CENTER,
                            "M", 48, COLOR_WHITE, NVG_ALIGN_CENTER);
}

void ui_draw_right_hud_infobox2(UIState *s) {
}

void ui_draw_right_hud(UIState *s) {
  ui_draw_right_hud_infobox1(s);
  ui_draw_right_hud_infobox2(s);
}

void ui_draw_infotext(UIState *s) {
#if 0
  int sidebar_fit_x = 0;
  char value[256];
  Params params;
  bool IsVagFulltimeLkaEnabled = params.getBool("IsVagFulltimeLkaEnabled");

  //Fit sidebar screen
  sidebar_fit_x = s->viz_rect.x + hud_left_2_x;

  snprintf(value, sizeof(value), "IsVagFulltimeLkaEnabled=%d", IsVagFulltimeLkaEnabled);
  ui_draw_hud_text(s, sidebar_fit_x, 5, value, 64, COLOR_YELLOW);
  snprintf(value, sizeof(value), "CoS enabled=%d, CoS.active=%d state=%d", s->scene.controls_state.getEnabled(), s->scene.controls_state.getActive(), (unsigned int)s->scene.controls_state.getState());
  ui_draw_hud_text(s, sidebar_fit_x, 700, value, 80, COLOR_YELLOW);
  snprintf(value, sizeof(value), "CaC.enabled=%d, CaC.active=%d CaC.availableFulltimeLka=%d", s->scene.car_control.getEnabled(), s->scene.car_control.getActive(), s->scene.car_control.getAvailableFulltimeLka());
  ui_draw_hud_text(s, sidebar_fit_x, 800, value, 80, COLOR_YELLOW);
  snprintf(value, sizeof(value), "CaS.CrS.enabled=%d, CaS.CrS.available=%d", s->scene.car_state.getCruiseState().getEnabled(), s->scene.car_state.getCruiseState().getAvailable());
  ui_draw_hud_text(s, sidebar_fit_x, 900, value, 80, COLOR_YELLOW);
#endif
}

void ui_draw_infobar(UIState *s) {
  const int x = s->viz_rect.x;
  const int y = s->viz_rect.bottom() - hud_infobar_h;
  const int w = s->viz_rect.w;
  const int text_x = w / 2 + x;
  const int text_y = y + 60;
  const bool brakeLights = s->scene.car_state.getBrakeLights();

  char infobar[100];
  // create time string
  char date_time[20];
  time_t rawtime = time(NULL);
  struct tm timeinfo;
  localtime_r(&rawtime, &timeinfo);
  strftime(date_time, sizeof(date_time),"%D %T", &timeinfo);
  snprintf(infobar, sizeof(infobar), "%s", date_time);

  nvgBeginPath(s->vg);
  nvgRect(s->vg, x, y, w, hud_infobar_h);
  nvgFillColor(s->vg, (brakeLights? COLOR_RED_ALPHA(200) : COLOR_BLACK_ALPHA(150)));
  nvgFill(s->vg);

  nvgFontSize(s->vg, 50);
  nvgFillColor(s->vg, COLOR_WHITE_ALPHA(200));
  nvgTextAlign(s->vg, NVG_ALIGN_CENTER);
  nvgText(s->vg, text_x, text_y, infobar, NULL);
}

void ui_draw_blinker(UIState *s) {
  // dp blinker, from kegman
  const int viz_blinker_w = 280;
  const int viz_blinker_x = s->viz_rect.centerX() - viz_blinker_w/2;
  const bool leftBlinker = s->scene.car_state.getLeftBlinker();
  const bool rightBlinker = s->scene.car_state.getRightBlinker();
  if(leftBlinker) {
    nvgBeginPath(s->vg);
    nvgMoveTo(s->vg, viz_blinker_x, s->viz_rect.y + header_h/4);
    nvgLineTo(s->vg, viz_blinker_x - viz_blinker_w/2, s->viz_rect.y + header_h/4 + header_h/4);
    nvgLineTo(s->vg, viz_blinker_x, s->viz_rect.y + header_h/2 + header_h/4);
    nvgClosePath(s->vg);
    nvgFillColor(s->vg, nvgRGBA(0,255,0,s->scene.ui_extend.blinker_blinkingrate>=60?190:30));
    nvgFill(s->vg);
  }
  if(rightBlinker) {
    nvgBeginPath(s->vg);
    nvgMoveTo(s->vg, viz_blinker_x+viz_blinker_w, s->viz_rect.y + header_h/4);
    nvgLineTo(s->vg, viz_blinker_x+viz_blinker_w + viz_blinker_w/2, s->viz_rect.y + header_h/4 + header_h/4);
    nvgLineTo(s->vg, viz_blinker_x+viz_blinker_w, s->viz_rect.y + header_h/2 + header_h/4);
    nvgClosePath(s->vg);
    nvgFillColor(s->vg, nvgRGBA(0,255,0,s->scene.ui_extend.blinker_blinkingrate>=60?190:30));
    nvgFill(s->vg);
  }
  if(leftBlinker || rightBlinker) {
    s->scene.ui_extend.blinker_blinkingrate -= 3;
    if(s->scene.ui_extend.blinker_blinkingrate<0) s->scene.ui_extend.blinker_blinkingrate = 120;
  }
}

void ui_draw_blindspot(UIState *s) {
  const int y = s->viz_rect.bottom() - hud_blindspot_w;
  const bool leftBlindspot = s->scene.car_state.getLeftBlindspot();
  const bool rightBlindspot = s->scene.car_state.getRightBlindspot();

  if (leftBlindspot) {
    const int left_x = s->viz_rect.x;
    nvgBeginPath(s->vg);
    nvgMoveTo(s->vg, left_x, y);
    nvgLineTo(s->vg, left_x, y+hud_blindspot_w);
    nvgLineTo(s->vg, left_x+hud_blindspot_w, y+hud_blindspot_w);
    nvgClosePath(s->vg);
    nvgFillColor(s->vg, COLOR_ORANGE_APPHA(200));
    nvgFill(s->vg);
  }

  if (rightBlindspot) {
    const int right_x = s->viz_rect.right();
    nvgBeginPath(s->vg);
    nvgMoveTo(s->vg, right_x, y);
    nvgLineTo(s->vg, right_x, y+hud_blindspot_w);
    nvgLineTo(s->vg, right_x-hud_blindspot_w, y+hud_blindspot_w);
    nvgClosePath(s->vg);
    nvgFillColor(s->vg, COLOR_ORANGE_APPHA(200));
    nvgFill(s->vg);
  }
}

static bool PreviousSpeedCameraDetected = false;
void ui_draw_speedcamera(UIState *s) {
  char speedLimit[16];
  char distance[16];

#if 0
  char value[64];
  int sidebar_fit_x = 0;
  //Fit sidebar screen
  sidebar_fit_x = s->viz_rect.x + hud_left_2_x;
  snprintf(value, sizeof(value), "Distance=%f, Angle=%f, VDirect=%d, SDirect=%d", \
           s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleDistance(), \
           s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleTrackAngle(), \
           (int) s->scene.speed_camera.getVehicleDirect(), \
           (int) s->scene.speed_camera.getSpeedCameraMapPosition().getDirect());
  ui_draw_hud_text(s, sidebar_fit_x, 5, value, 64, COLOR_YELLOW);

  snprintf(value, sizeof(value), "V latitude=%f, longitude=%f", \
           s->scene.speed_camera.getVehicleLatitude(), \
           s->scene.speed_camera.getVehicleLongitude());
  ui_draw_hud_text(s, sidebar_fit_x, 800, value, 80, COLOR_YELLOW);

  snprintf(value, sizeof(value), "C latitude=%f, longitude=%f", \
           s->scene.speed_camera.getSpeedCameraMapPosition().getLatitude(), \
           s->scene.speed_camera.getSpeedCameraMapPosition().getLongitude());
  ui_draw_hud_text(s, sidebar_fit_x, 900, value, 80, COLOR_YELLOW);
#endif

  if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 25.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_25");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 30.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_30");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 40.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_40");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 50.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_50");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 60.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_60");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 70.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_70");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 80.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_80");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 90.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_90");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 100.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_100");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getSpeedLimitation() == 110.0) {
    snprintf(speedLimit, sizeof(speedLimit), "speed_limit_110");
  }

  if(s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleDistance() < 0.2) {
    snprintf(distance, sizeof(distance), "100M");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleDistance() < 0.3) {
    snprintf(distance, sizeof(distance), "200M");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleDistance() < 0.4) {
    snprintf(distance, sizeof(distance), "300M");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleDistance() < 0.5) {
    snprintf(distance, sizeof(distance), "400M");
  } else if(s->scene.speed_camera.getSpeedCameraMapPosition().getVehicleDistance() < 0.6) {
    snprintf(distance, sizeof(distance), "500M");
  }

  if(s->scene.speed_camera.getSpeedCameraDetected()) {
    ui_draw_speed_image(s, 1650, 500, 200, 200, speedLimit);
    ui_draw_hud_text(s, 1700, 700, distance, 64, COLOR_YELLOW);
    if(PreviousSpeedCameraDetected == false) {
      //PONTEST s->sound->play(cereal::CarControl::HUDControl::AudibleAlert(4));
    }
  }
  PreviousSpeedCameraDetected = s->scene.speed_camera.getSpeedCameraDetected();
}

//===== draw hud =====
void ui_draw_hud(UIState *s) {
  bool IsVagInfoboxEnabled = true;
  bool IsVagInfobarEnabled = true;
  bool IsVagBlinkerEnabled = true;
  bool IsVagBlindspotEnabled = true;
  bool IsVagSpeedCameraEnabled = true;
  bool IsVagDevelopModeEnabled = true;
  Params params;

  IsVagInfoboxEnabled = params.getBool("IsVagInfoboxEnabled");
  IsVagInfobarEnabled = params.getBool("IsVagInfobarEnabled");
  IsVagBlinkerEnabled = params.getBool("IsVagBlinkerEnabled");
  IsVagBlindspotEnabled = params.getBool("IsVagBlindspotEnabled");
  IsVagSpeedCameraEnabled = params.getBool("IsVagSpeedCameraEnabled");
  IsVagDevelopModeEnabled = params.getBool("IsVagDevelopModeEnabled");

  if(IsVagInfoboxEnabled) {
    //ui_draw_top_hud(s);
    //ui_draw_bottom_hud(s);
    ui_draw_left_hud(s);
    //ui_draw_right_hud(s);
  }
  if(IsVagDevelopModeEnabled) {
   ui_draw_infotext(s);
  }
  if(IsVagInfobarEnabled) {
    ui_draw_infobar(s);
  }
  if(IsVagBlindspotEnabled) {
    ui_draw_blindspot(s);
  }
  if(IsVagBlinkerEnabled) {
    ui_draw_blinker(s);
  }
  if(IsVagSpeedCameraEnabled) {
    ui_draw_speedcamera(s);
  }
}
