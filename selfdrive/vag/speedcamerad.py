#!/usr/bin/env python3
import os
import math
import time
import csv
from cereal import car, log

from common.numpy_fast import clip
from common.realtime import sec_since_boot, config_realtime_process, Priority, Ratekeeper, DT_CTRL
from common.profiler import Profiler
from common.params import Params, put_nonblocking
import cereal.messaging as messaging
from selfdrive.config import Conversions as CV
from selfdrive.swaglog import cloudlog
from selfdrive.boardd.boardd import can_list_to_can_capnp
from selfdrive.car.car_helpers import get_car, get_startup_event, get_one_can
from selfdrive.controls.lib.lane_planner import CAMERA_OFFSET
from selfdrive.controls.lib.drive_helpers import update_v_cruise, initialize_v_cruise
from selfdrive.controls.lib.longcontrol import LongControl, STARTING_TARGET_SPEED
from selfdrive.controls.lib.latcontrol_pid import LatControlPID
from selfdrive.controls.lib.latcontrol_indi import LatControlINDI
from selfdrive.controls.lib.latcontrol_lqr import LatControlLQR
from selfdrive.controls.lib.latcontrol_angle import LatControlAngle
from selfdrive.controls.lib.events import Events, ET
from selfdrive.controls.lib.alertmanager import AlertManager
from selfdrive.controls.lib.vehicle_model import VehicleModel
from selfdrive.controls.lib.longitudinal_planner import LON_MPC_STEP
from selfdrive.locationd.calibrationd import Calibration
from selfdrive.hardware import HARDWARE, TICI

LDW_MIN_SPEED = 31 * CV.MPH_TO_MS
LANE_DEPARTURE_THRESHOLD = 0.1
STEER_ANGLE_SATURATION_TIMEOUT = 1.0 / DT_CTRL
STEER_ANGLE_SATURATION_THRESHOLD = 2.5  # Degrees

SIMULATION = "SIMULATION" in os.environ
NOSENSOR = "NOSENSOR" in os.environ
IGNORE_PROCESSES = set(["rtshield", "uploader", "deleter", "loggerd", "logmessaged", "tombstoned", "logcatd", "proclogd", "clocksd", "updated", "timezoned", "manage_athenad"])

ThermalStatus = log.DeviceState.ThermalStatus
State = log.ControlsState.OpenpilotState
PandaType = log.PandaState.PandaType
Desire = log.LateralPlan.Desire
LaneChangeState = log.LateralPlan.LaneChangeState
LaneChangeDirection = log.LateralPlan.LaneChangeDirection
EventName = car.CarEvent.EventName
GearShifter = car.CarState.GearShifter

class MapPoint:
  def __init__(self, Longitude, Latitude):
    self.Longitude = Longitude
    self.Latitude = Latitude

class CameraMapPoint(MapPoint):
  def __init__(self, Longitude, Latitude, Direct, SpeedLimit):
    super().__init__(Longitude, Latitude)
    self.Direct = Direct
    self.SpeedLimit = SpeedLimit

  def distance(self, MapPoint):
    return MapPoint.Longitude - self.Longitude, MapPoint.Latitude - self.Latitude

class SpeedCamera:
  def __init__(self):
    print("[PONTEST][speedcamerad.py][__init__()]")
    #self.gps = messaging.sub_sock('gpsLocationExternal')
    #self.SpeedCameraMapPointList = []
    #self.PassiveLatitudeList = []
    #self.NagtiveLatitudeList = []
    #self.PassiveLongitudeList = []
    #self.NagtiveLongitudeList = []
    path = '/data/openpilot/selfdrive/vag/speedcamera_csv/NPA_TD1.csv'
    try:
      f = open(path, 'r')
      rows = csv.reader(f, delimiter=',')
      for row in rows:
        print(row)
    except:
      print('ERROR: can not found ' + path)
      exit(1)


  def update_events(self):
    print("[PONTEST][speedcamerad.py][update_events()]")
    #print("[PONTEST][speedcamerad.py][update_events() latitude=]", self.gps.latitude, " longitude=", self.gps.longitude)

  def speedcamerad_thread(self):
    print("[PONTEST][speedcamerad.py][speedcamerad_thread()] self.gps")
    while True:
      self.update_events()
      time.sleep(1)

def main():
  print("[PONTEST][speedcamerad.py][main()]")
  speedcamera = SpeedCamera()
  speedcamera.speedcamerad_thread()


if __name__ == "__main__":
  main()
