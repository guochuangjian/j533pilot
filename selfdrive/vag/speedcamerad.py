#!/usr/bin/env python3
import os
import math
import time
import csv
from decimal import Decimal
from cereal import car, log
import cereal.messaging as messaging

from common.numpy_fast import clip
from common.realtime import sec_since_boot, config_realtime_process, Priority, Ratekeeper, DT_CTRL
from common.profiler import Profiler
from common.params import Params, put_nonblocking
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

class MapPosition:
  def __init__(self, Latitude, Longitude, Direct):
    self.Latitude = Latitude
    self.Longitude = Longitude
    self.Direct = Direct

class SpeedCameraMapPosition(MapPosition):
  def __init__(self, Latitude, Longitude, Direct, SpeedLimit, RoadType, Distance=-1, Angle=-1):
    super().__init__(Latitude, Longitude, Direct)
    self.SpeedLimit = SpeedLimit
    self.RoadType = RoadType
    self.Distance = Distance
    self.Angle = Angle

class SpeedCamera:
  def __init__(self):
    print("[PONTEST][speedcamerad.py][__init__()]")
    #self.gps = messaging.sub_sock('gpsLocationExternal')
    self.VehicleMapPositionList = []
    self.SpeedCameraMapPositionList = []
    self.ConcentricLayer1PositionList = [] #under 5km
    self.ConcentricLayer2PositionList = [] #5~10km
    self.ConcentricLayer3PositionList = [] #over 10km
    self.TestItemIndex = 0

    #Add speed camera map position
    SpeedCamPath = '/data/openpilot/selfdrive/vag/speedcamera_csv/TaiwanSpeedCamera.csv'
    try:
      f = open(SpeedCamPath, 'r')
      rows = csv.reader(f, delimiter=',')
      #print(time.ctime())
      for row in rows:
        self.SpeedCameraMapPositionList.append(SpeedCameraMapPosition(row[0], row[1], row[2], row[3], row[4]))
        #print(row[0], row[1], row[2], row[3], row[4])
        #print("SpeedCameraMapPositionList=", self.SpeedCameraMapPositionList[-1].Latitude, \
        #                                  self.SpeedCameraMapPositionList[-1].Longitude, \
        #                                  self.SpeedCameraMapPositionList[-1].Direct, \
        #                                  self.SpeedCameraMapPositionList[-1].SpeedLimit, \
        #                                  self.SpeedCameraMapPositionList[-1].RoadType)
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
        if not self.VehicleMapPositionList:
          #print("VehicleMapPositionList<=0")
          self.VehicleMapPositionList.append(MapPosition(row[0], row[1], ''))
        else:
          #print("VehicleMapPositionList>0")
          #print(row[0], self.VehicleMapPositionList[-1].Latitude, row[1], self.VehicleMapPositionList[-1].Longitude)
          LatitudeDistance = Decimal(row[0]) - Decimal(self.VehicleMapPositionList[-1].Latitude)
          LongitudeDistance = Decimal(row[1]) - Decimal(self.VehicleMapPositionList[-1].Longitude)

          if abs(LongitudeDistance) > abs(LatitudeDistance):
            if LongitudeDistance > 0:
              #print("E")
              self.VehicleMapPositionList.append(MapPosition(row[0], row[1], 'E'))
            else:
              #print("W")
              self.VehicleMapPositionList.append(MapPosition(row[0], row[1], 'W'))

          else:
            if LatitudeDistance > 0:
              #print("N")
              self.VehicleMapPositionList.append(MapPosition(row[0], row[1], 'N'))
            else:
              #print("S")
              self.VehicleMapPositionList.append(MapPosition(row[0], row[1], 'S'))

        #print(row[0], row[1])
        #print("VehicleMapPositionList=", self.VehicleMapPositionList[-1].Latitude, self.VehicleMapPositionList[-1].Longitude, self.VehicleMapPositionList[-1].Direct)

      fd.close()
      print(time.ctime())
    except:
      print('ERROR: can not found ' + VehiclePath)
      exit(1)

    #Sort SpeedCameraMapPositionList
    self.SpeedCameraMapPositionList.sort(key = lambda s: s.Longitude)
    self.SpeedCameraMapPositionList.sort(key = lambda s: s.Latitude, reverse = True)

    #for item in self.SpeedCameraMapPositionList:
    #  print("item=", item.Latitude, \
    #                 item.Longitude, \
    #                 item.Direct, \
    #                 item.SpeedLimit, \
    #                 item.RoadType)

    #Divide concentric layer
    #print("VehiclePosition=", self.VehicleMapPositionList[0].Latitude, self.VehicleMapPositionList[0].Longitude)
    try:
      for MapPositionItem in self.SpeedCameraMapPositionList:
        #print("test", self.VehicleMapPositionList[0].Latitude, \
        #              self.VehicleMapPositionList[0].Longitude, \
        #              MapPositionItem.Latitude, \
        #              MapPositionItem.Longitude)
        distance1 = self.calculate_gps_radian_great_circle_distance(self.VehicleMapPositionList[0].Latitude, \
                                                                    self.VehicleMapPositionList[0].Longitude, \
                                                                    MapPositionItem.Latitude, \
                                                                    MapPositionItem.Longitude)
        distance2 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[0].Latitude, \
                                                                         self.VehicleMapPositionList[0].Longitude, \
                                                                         MapPositionItem.Latitude, \
                                                                         MapPositionItem.Longitude)
        #print("SpeedCamPosition=", MapPositionItem.Latitude, MapPositionItem.Longitude, distance1, distance2)
    except:
      print("Divide concentric layer fail")
      exit(1)


  def calculate_gps_radian_great_circle_distance(self, Latitude1, Longitude1, Latitude2, Longitude2):
    #print("[PONTEST][speedcamerad.py][calculate_gps_radian_great_circle_distance()]")
    R = 6371 # Earth's radius (km)
    d = R * math.acos(math.sin(Decimal(Latitude1)) * math.sin(Decimal(Latitude2)) + \
          math.cos(Decimal(Latitude1)) * math.cos(Decimal(Latitude2)) * math.cos(Decimal(Longitude2) - Decimal(Longitude1)))
    return d

  def calculate_gps_radian_haversine_formula_distance(self, Latitude1, Longitude1, Latitude2, Longitude2):
    #print("[PONTEST][speedcamerad.py][calculate_gps_radian_haversine_formula_distance()]")
    Latitude1Radian = Decimal(Latitude1) * Decimal(math.pi) / 180
    Longitude1Radian = Decimal(Longitude1) * Decimal(math.pi) / 180
    Latitude2Radian = Decimal(Latitude2) * Decimal(math.pi) / 180
    Longitude2Radian = Decimal(Longitude2) * Decimal(math.pi) / 180
    #print("[PONTEST][speedcamerad.py][calculate_gps_radian_haversine_formula_distance()] 2")

    a = Decimal(Latitude1Radian) - Decimal(Latitude2Radian)
    b = Decimal(Longitude1Radian) - Decimal(Longitude2Radian)
    #cal = math.sqrt(math.pow(Decimal(a), 2))
    cal = 2 * math.asin(math.sqrt(math.pow(math.sin(Decimal(a) / 2), 2) + \
          math.cos(Decimal(Latitude1Radian)) * math.cos(Decimal(Latitude2Radian)) * math.pow(math.sin(Decimal(b) / 2), 2))) * 6378.137
    #distance = math.round(Decimal(cal) * 10000, 2) / 10000
    #distance = Decimal(cal) * 10000 / 10000
    #print("[PONTEST][speedcamerad.py][calculate_gps_radian_haversine_formula_distance()] 3")
    return cal

  def calculate_position_angle(self, Distance1, Distance2, Distance3):
    #print("[PONTEST][speedcamerad.py][calculate_position_angle()]")
    #angle = (math.pow(Decimal(Distance2), 2) + math.pow(Decimal(Distance3), 2) - math.pow(Decimal(Distance1), 2)) / (2 * Decimal(Distance2) * Decimal(Distance3))
    if Distance2 > 0 and Distance3 > 0:
      angle = math.degrees(math.acos((math.pow(Decimal(Distance2), 2) + math.pow(Decimal(Distance3), 2) - math.pow(Decimal(Distance1), 2)) / (2 * (Distance2) * (Distance3))))
    else:
      angle = 361
    return angle
    

  def recalculate_concentric_layer_all(self):
    print("[PONTEST][speedcamerad.py][recalculate_concentric_layer_all()]")
    try:
      for MapPositionItem in self.SpeedCameraMapPositionList:
        #print("test", self.VehicleMapPositionList[self.TestItemIndex].Latitude, \
        #              self.VehicleMapPositionList[self.TestItemIndex].Longitude, \
        #              MapPositionItem.Latitude, \
        #              MapPositionItem.Longitude)
        distance1 = self.calculate_gps_radian_great_circle_distance(self.VehicleMapPositionList[self.TestItemIndex].Latitude, \
                                                                    self.VehicleMapPositionList[self.TestItemIndex].Longitude, \
                                                                    MapPositionItem.Latitude, \
                                                                    MapPositionItem.Longitude)
        distance2 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[self.TestItemIndex].Latitude, \
                                                                         self.VehicleMapPositionList[self.TestItemIndex].Longitude, \
                                                                         MapPositionItem.Latitude, \
                                                                         MapPositionItem.Longitude)
        
        if self.TestItemIndex > 0:
          distance3 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[self.TestItemIndex - 1].Latitude, \
                                                                           self.VehicleMapPositionList[self.TestItemIndex - 1].Longitude, \
                                                                           MapPositionItem.Latitude, \
                                                                           MapPositionItem.Longitude)
          distance4 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[self.TestItemIndex - 1].Latitude, \
                                                                           self.VehicleMapPositionList[self.TestItemIndex - 1].Longitude, \
                                                                           self.VehicleMapPositionList[self.TestItemIndex].Latitude, \
                                                                           self.VehicleMapPositionList[self.TestItemIndex].Longitude)
          angle = self.calculate_position_angle(distance2, distance3, distance4)
        else:
          angle = -1

        #print("SpeedCamPosition=", self.VehicleMapPositionList[self.TestItemIndex].Latitude, self.VehicleMapPositionList[self.TestItemIndex].Longitude, MapPositionItem.Latitude, MapPositionItem.Longitude, distance1, distance2, angle)
        if distance2 <= 5:
          self.ConcentricLayer1PositionList.append(SpeedCameraMapPosition(MapPositionItem.Latitude, MapPositionItem.Longitude, MapPositionItem.Direct, MapPositionItem.SpeedLimit, MapPositionItem.RoadType, distance2, angle))
        elif distance2 > 5 and distance2 <=10:
          self.ConcentricLayer2PositionList.append(SpeedCameraMapPosition(MapPositionItem.Latitude, MapPositionItem.Longitude, MapPositionItem.Direct, MapPositionItem.SpeedLimit, MapPositionItem.RoadType, distance2, angle))
        else:
          self.ConcentricLayer3PositionList.append(SpeedCameraMapPosition(MapPositionItem.Latitude, MapPositionItem.Longitude, MapPositionItem.Direct, MapPositionItem.SpeedLimit, MapPositionItem.RoadType, distance2, angle))

      self.ConcentricLayer1PositionList.sort(key = lambda s: s.Distance)
      self.ConcentricLayer2PositionList.sort(key = lambda s: s.Distance)
      self.ConcentricLayer3PositionList.sort(key = lambda s: s.Distance)
    except:
      print("recalculate_concentric_layer_all fail")
      exit(1)
    print("len(self.ConcentricLayer1PositionList)=", len(self.ConcentricLayer1PositionList))
    for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
      print("ConcentricLayer1Item=", ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)

    print("len(self.ConcentricLayer2PositionList)=", len(self.ConcentricLayer2PositionList))
    #for ConcentricLayer2Item in self.ConcentricLayer2PositionList:
    #  print("ConcentricLayer2Item=", ConcentricLayer2Item.Latitude, ConcentricLayer2Item.Longitude, ConcentricLayer2Item.Direct, ConcentricLayer2Item.SpeedLimit, ConcentricLayer2Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)

    print("len(self.ConcentricLayer3PositionList)=", len(self.ConcentricLayer3PositionList))
    #for ConcentricLayer3Item in self.ConcentricLayer3PositionList:
    #  print("ConcentricLayer3Item=", ConcentricLayer3Item.Latitude, ConcentricLayer3Item.Longitude, ConcentricLayer3Item.Direct, ConcentricLayer3Item.SpeedLimit, ConcentricLayer3Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)


  def recalculate_concentric_layer1(self):
    print("[PONTEST][speedcamerad.py][recalculate_concentric_layer1()]")
    for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
      distance2 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[self.TestItemIndex].Latitude, \
                                                                       self.VehicleMapPositionList[self.TestItemIndex].Longitude, \
                                                                       ConcentricLayer1Item.Latitude, \
                                                                       ConcentricLayer1Item.Longitude)
      if self.TestItemIndex > 0:
        distance3 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[self.TestItemIndex - 1].Latitude, \
                                                                         self.VehicleMapPositionList[self.TestItemIndex - 1].Longitude, \
                                                                         ConcentricLayer1Item.Latitude, \
                                                                         ConcentricLayer1Item.Longitude)
        distance4 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[self.TestItemIndex - 1].Latitude, \
                                                                         self.VehicleMapPositionList[self.TestItemIndex - 1].Longitude, \
                                                                         self.VehicleMapPositionList[self.TestItemIndex].Latitude, \
                                                                         self.VehicleMapPositionList[self.TestItemIndex].Longitude)
        angle = self.calculate_position_angle(distance2, distance3, distance4)
      else:
        angle = -1

      if distance2 > 5:
        self.ConcentricLayer1PositionList.remove(ConcentricLayer1Item)
      else:
        ConcentricLayer1Item.Distance = distance2
        ConcentricLayer1Item.Angle = angle

    self.ConcentricLayer1PositionList.sort(key = lambda s: s.Distance)
    print("len(self.ConcentricLayer1PositionList)=", len(self.ConcentricLayer1PositionList))
    for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
      print("ConcentricLayer1Item=", ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)

  def vehicle_track_simulate(self):
    print("[PONTEST][speedcamerad.py][vehicle_track_simulate()]")
    if self.TestItemIndex < 1900:
      print("[PONTEST][speedcamerad.py][vehicle_track_simulate()] simulate start")
      self.TestItemIndex = self.TestItemIndex + 1

      if len(self.ConcentricLayer1PositionList) < 5:
        self.recalculate_concentric_layer_all()
      else:
        self.recalculate_concentric_layer1()

      #for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
      #  print("ConcentricLayer1Item=", self.TestItemIndex, self.VehicleMapPositionList[self.TestItemIndex].Latitude, self.VehicleMapPositionList[self.TestItemIndex].Longitude, ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance)
      print("SpeedCamItem=", self.TestItemIndex, self.VehicleMapPositionList[self.TestItemIndex].Latitude, self.VehicleMapPositionList[self.TestItemIndex].Longitude, self.ConcentricLayer1PositionList[0].Latitude, self.ConcentricLayer1PositionList[0].Longitude, self.ConcentricLayer1PositionList[0].Direct, self.ConcentricLayer1PositionList[0].SpeedLimit, self.ConcentricLayer1PositionList[0].RoadType, self.ConcentricLayer1PositionList[0].Distance, self.ConcentricLayer1PositionList[0].Angle)
      print("[PONTEST][speedcamerad.py][vehicle_track_simulate()] simulate end")


  def update_events(self, VehicleLatitude, VehicleLongitude, VehicleSpeed):
    print("[PONTEST][speedcamerad.py][update_events()]", VehicleLatitude, VehicleLongitude, VehicleSpeed)

  def speedcamerad_thread(self):
    print("[PONTEST][speedcamerad.py][speedcamerad_thread()] self.gps")
    location_sock = messaging.sub_sock('gpsLocationExternal')

    while True:
      location = messaging.recv_sock(location_sock)
      location = (strip_deprecated_keys(location.gpsLocationExternal.to_dict()) if location else None)
      if location:
        print("location=", location.latitude, location.longitude, location.altitude, location.speed, location.timestamp, location.source)
        self.update_events(location.latitude, location.longitude, location.speed)
      time.sleep(100)

def main():
  print("[PONTEST][speedcamerad.py][main()]")
  speedcamera = SpeedCamera()
  speedcamera.vehicle_track_simulate()
  speedcamera.speedcamerad_thread()


if __name__ == "__main__":
  main()
