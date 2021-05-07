#!/usr/bin/env python3
import os
import math
import time
import csv
from decimal import Decimal
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
  def __init__(self, Latitude, Longitude, Direct):
    self.Latitude = Latitude
    self.Longitude = Longitude
    self.Direct = Direct

class SpeedCameraMapPoint(MapPoint):
  def __init__(self, Latitude, Longitude, Direct, SpeedLimit, RoadType):
    super().__init__(Latitude, Longitude, Direct)
    self.SpeedLimit = SpeedLimit
    self.RoadType = RoadType

class SpeedCamera:
  def __init__(self):
    print("[PONTEST][speedcamerad.py][__init__()]")
    #self.gps = messaging.sub_sock('gpsLocationExternal')
    self.VehicleMapPointList = []
    self.SpeedCameraMapPointList = []
    self.ConcentricLayer1PointList = [] #under 5km
    self.ConcentricLayer2PointList = [] #5~10km
    self.ConcentricLayer3PointList = [] #over 10km

    #Add speed camera map point
    SpeedCamPath = '/data/openpilot/selfdrive/vag/speedcamera_csv/TaiwanSpeedCamera.csv'
    try:
      f = open(SpeedCamPath, 'r')
      rows = csv.reader(f, delimiter=',')
      #print(time.ctime())
      for row in rows:
        self.SpeedCameraMapPointList.append(SpeedCameraMapPoint(row[0], row[1], row[2], row[3], row[4]))
        #print(row[0], row[1], row[2], row[3], row[4])
        #print("SpeedCameraMapPointList=", self.SpeedCameraMapPointList[-1].Latitude, \
        #                                  self.SpeedCameraMapPointList[-1].Longitude, \
        #                                  self.SpeedCameraMapPointList[-1].Direct, \
        #                                  self.SpeedCameraMapPointList[-1].SpeedLimit, \
        #                                  self.SpeedCameraMapPointList[-1].RoadType)
      f.close()
      #print(time.ctime())
    except:
      print('ERROR: can not found ' + SpeedCamPath)
      exit(1)

    #Add vehicle tracke data
    VehiclePath = '/data/openpilot/selfdrive/vag/speedcamera_csv/GpsTrack.csv'
    try:
      fd = open(VehiclePath, 'r')
      rows = csv.reader(fd, delimiter=',')
      #print(time.ctime())
      for row in rows:
        if not self.VehicleMapPointList:
          #print("VehicleMapPointList<=0")
          self.VehicleMapPointList.append(MapPoint(row[0], row[1], ''))
        else:
          #print("VehicleMapPointList>0")
          #print(row[0], self.VehicleMapPointList[-1].Latitude, row[1], self.VehicleMapPointList[-1].Longitude)
          LatitudeDistance = Decimal(row[0]) - Decimal(self.VehicleMapPointList[-1].Latitude)
          LongitudeDistance = Decimal(row[1]) - Decimal(self.VehicleMapPointList[-1].Longitude)

          if abs(LongitudeDistance) > abs(LatitudeDistance):
            if LongitudeDistance > 0:
              #print("E")
              self.VehicleMapPointList.append(MapPoint(row[0], row[1], 'E'))
            else:
              #print("W")
              self.VehicleMapPointList.append(MapPoint(row[0], row[1], 'W'))

          else:
            if LatitudeDistance > 0:
              #print("N")
              self.VehicleMapPointList.append(MapPoint(row[0], row[1], 'N'))
            else:
              #print("S")
              self.VehicleMapPointList.append(MapPoint(row[0], row[1], 'S'))

        #print(row[0], row[1])
        print("VehicleMapPointList=", self.VehicleMapPointList[-1].Latitude, self.VehicleMapPointList[-1].Longitude, self.VehicleMapPointList[-1].Direct)

      fd.close()
      print(time.ctime())
    except:
      print('ERROR: can not found ' + VehiclePath)
      exit(1)

    #Sort SpeedCameraMapPointList
    self.SpeedCameraMapPointList.sort(key = lambda s: s.Longitude)
    self.SpeedCameraMapPointList.sort(key = lambda s: s.Latitude, reverse = True)

    for item in self.SpeedCameraMapPointList:
      print("item=", item.Latitude, \
                     item.Longitude, \
                     item.Direct, \
                     item.SpeedLimit, \
                     item.RoadType)

    #Divide concentric layer
    print("VehiclePosition=", self.VehicleMapPointList[0].Latitude, self.VehicleMapPointList[0].Longitude)
    try:
      for MapPointItem in self.SpeedCameraMapPointList:
        print("test", self.VehicleMapPointList[0].Latitude, \
                      self.VehicleMapPointList[0].Longitude, \
                      MapPointItem.Latitude, \
                      MapPointItem.Longitude)
        distance = self.get_gps_radian_distance(self.VehicleMapPointList[0].Latitude, \
                                                self.VehicleMapPointList[0].Longitude, \
                                                MapPointItem.Latitude, \
                                                MapPointItem.Longitude)
        print("SpeedCamPosition=", MapPointItem.Latitude, MapPointItem.Longitude, distance)
    except:
      print("Divide concentric layer fail")
      exit(1)


  def get_gps_radian_distance(self, Latitude1, Longitude1, Latitude2, Longitude2):
    print("[PONTEST][speedcamerad.py][get_gps_radian_distance()]")
    R = 6371 # Earth's radius (km)
    d = R * math.acos(math.sin(Decimal(Latitude1)) * math.sin(Decimal(Latitude2)) + \
          math.cos(Decimal(Latitude1)) * math.cos(Decimal(Latitude2)) * math.cos(Decimal(Longitude2) - Decimal(Longitude1)))
    return d

  def get_gps_degree_distance(self):
    print("[PONTEST][speedcamerad.py][get_gps_degree_distance()]")

  def update_events(self):
    print("[PONTEST][speedcamerad.py][update_events()]")

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
