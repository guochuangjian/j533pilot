#!/usr/bin/env python3
import os
import math
import time
import csv
import sys
from decimal import Decimal
from cereal import car, log
import cereal.messaging as messaging

SpeedDirect = log.SpeedCamera.SpeedDirect
RoadType = log.SpeedCamera.SpeedCameraMapPosition.RoadType

class MapPosition:
  def __init__(self, Latitude, Longitude, Direct):
    self.Latitude = Latitude
    self.Longitude = Longitude
    self.Direct = Direct


class SpeedCameraMapPosition(MapPosition):
  def __init__(self, Latitude, Longitude, Direct, SpeedLimit, Type, Distance=-1, Angle=-1):
    super().__init__(Latitude, Longitude, Direct)
    self.SpeedLimit = SpeedLimit
    self.RoadType = Type
    self.Distance = Distance
    self.Angle = Angle


class SpeedCamera:
  def __init__(self):
    print("[PONTEST][speedcamerad.py][__init__()]")
    self.VehiclePreviousLatitude = (360.1)
    self.VehiclePreviousLongitude = (360.1)
    self.VehiclePreviousSpeedCameraDistance = 0;
    self.VehicleToCameraDirect = SpeedDirect.u
    self.SpeedCameraDetected = False
    self.VehicleMapPositionList = []
    self.SpeedCameraMapPositionList = []
    self.ConcentricLayer1PositionList = [] #under 5km
    self.ConcentricLayer2PositionList = [] #5~10km
    self.ConcentricLayer3PositionList = [] #over 10km
    self.sm = messaging.SubMaster(['gpsLocationExternal'])
    self.pm = messaging.PubMaster(['speedCamera'])
    self.TestItemIndex = 0

    #Add speed camera map position
    SpeedCamPath = '/data/openpilot/selfdrive/vag/speedcamera_csv/TaiwanSpeedCamera.csv'
    f = open(SpeedCamPath, 'r')
    rows = csv.reader(f, delimiter=',')
    for row in rows:
      if (row[2] == 'N'):
        Direct = SpeedDirect.n
      elif (row[2] == 'S'):
        Direct = SpeedDirect.s
      elif (row[2] == 'E'):
        Direct = SpeedDirect.e
      elif (row[2] == 'W'):
        Direct = SpeedDirect.w
      elif (row[2] == 'D'):
        Direct = SpeedDirect.d
      else:
        Direct = SpeedDirect.d

      if (row[4] == 'road'):
        Type = RoadType.road
      elif (row[4] == 'freeway'):
        Type = RoadType.freeway
      elif (row[4] == 'highway'):
        Type = RoadType.highway
      else:
        Type = RoadType.road

      self.SpeedCameraMapPositionList.append(SpeedCameraMapPosition(row[0], row[1], Direct, row[3], Type))
      #self.SpeedCameraMapPositionList.append(SpeedCameraMapPosition(row[0], row[1], row[2], row[3], row[4]))
      #print(row[0], row[1], row[2], row[3], row[4])
      #print("SpeedCameraMapPositionList=", self.SpeedCameraMapPositionList[-1].Latitude, \
      #                                  self.SpeedCameraMapPositionList[-1].Longitude, \
      #                                  self.SpeedCameraMapPositionList[-1].Direct, \
      #                                  self.SpeedCameraMapPositionList[-1].SpeedLimit, \
      #                                  self.SpeedCameraMapPositionList[-1].RoadType)

    f.close()


    #Sort SpeedCameraMapPositionList
    self.SpeedCameraMapPositionList.sort(key = lambda s: s.Longitude)
    self.SpeedCameraMapPositionList.sort(key = lambda s: s.Latitude, reverse = True)

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


  def recalculate_concentric_layer_all(self, VehicleLatitude, VehicleLongitude, VehicleSpeed):
    #print("[PONTEST][speedcamerad.py][recalculate_concentric_layer_all()]")

    del self.ConcentricLayer1PositionList[:]
    del self.ConcentricLayer2PositionList[:]
    del self.ConcentricLayer3PositionList[:]
    #try:
    for MapPositionItem in self.SpeedCameraMapPositionList:
      if (float(MapPositionItem.Latitude) > float(VehicleLatitude) - 0.025) and \
         (float(MapPositionItem.Latitude) < float(VehicleLatitude) + 0.025) and \
         (float(MapPositionItem.Longitude) > float(VehicleLongitude) - 0.025) and \
         (float(MapPositionItem.Longitude) < float(VehicleLongitude) + 0.025):
        self.ConcentricLayer1PositionList.append(SpeedCameraMapPosition(MapPositionItem.Latitude, MapPositionItem.Longitude, MapPositionItem.Direct, MapPositionItem.SpeedLimit, MapPositionItem.RoadType, -1, -1))
      elif (float(MapPositionItem.Latitude) > float(VehicleLatitude) - 0.5) and \
           (float(MapPositionItem.Latitude) < float(VehicleLatitude) + 0.5) and \
           (float(MapPositionItem.Longitude) > float(VehicleLongitude) - 0.5) and \
           (float(MapPositionItem.Longitude) < float(VehicleLongitude) + 0.5):
        self.ConcentricLayer2PositionList.append(SpeedCameraMapPosition(MapPositionItem.Latitude, MapPositionItem.Longitude, MapPositionItem.Direct, MapPositionItem.SpeedLimit, MapPositionItem.RoadType, -1, -1))
      else:
        self.ConcentricLayer3PositionList.append(SpeedCameraMapPosition(MapPositionItem.Latitude, MapPositionItem.Longitude, MapPositionItem.Direct, MapPositionItem.SpeedLimit, MapPositionItem.RoadType, -1, -1))

    self.ConcentricLayer1PositionList.sort(key = lambda s: s.Distance)
    self.ConcentricLayer2PositionList.sort(key = lambda s: s.Distance)
    self.ConcentricLayer3PositionList.sort(key = lambda s: s.Distance)
    #except:
    #  print('[SpeedCamera::recalculate_concentric_layer_all()] Execption ')
    #  exit(1)
    print("[PONTEST][speedcamerad.py][recalculate_concentric_layer_all()] len(self.ConcentricLayer1PositionList)=", len(self.ConcentricLayer1PositionList))
    #for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
    #  print("ConcentricLayer1Item=", ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)

    print("[PONTEST][speedcamerad.py][recalculate_concentric_layer_all()] len(self.ConcentricLayer2PositionList)=", len(self.ConcentricLayer2PositionList))
    #for ConcentricLayer2Item in self.ConcentricLayer2PositionList:
    #  print("ConcentricLayer2Item=", ConcentricLayer2Item.Latitude, ConcentricLayer2Item.Longitude, ConcentricLayer2Item.Direct, ConcentricLayer2Item.SpeedLimit, ConcentricLayer2Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)

    print("[PONTEST][speedcamerad.py][recalculate_concentric_layer_all()] len(self.ConcentricLayer3PositionList)=", len(self.ConcentricLayer3PositionList))
    #for ConcentricLayer3Item in self.ConcentricLayer3PositionList:
    #  print("ConcentricLayer3Item=", ConcentricLayer3Item.Latitude, ConcentricLayer3Item.Longitude, ConcentricLayer3Item.Direct, ConcentricLayer3Item.SpeedLimit, ConcentricLayer3Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)


  def recalculate_concentric_layer1(self, VehicleLatitude, VehicleLongitude, VehicleSpeed):
    #print("[PONTEST][speedcamerad.py][recalculate_concentric_layer1()]")

    #del self.ConcentricLayer1PositionList[:]
    for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
      distance2 = self.calculate_gps_radian_haversine_formula_distance(VehicleLatitude, \
                                                                       VehicleLongitude, \
                                                                       ConcentricLayer1Item.Latitude, \
                                                                       ConcentricLayer1Item.Longitude)
      if self.TestItemIndex > 0:
        distance3 = self.calculate_gps_radian_haversine_formula_distance(self.VehiclePreviousLatitude, \
                                                                         self.VehiclePreviousLongitude, \
                                                                         ConcentricLayer1Item.Latitude, \
                                                                         ConcentricLayer1Item.Longitude)
        distance4 = self.calculate_gps_radian_haversine_formula_distance(self.VehiclePreviousLatitude, \
                                                                         self.VehiclePreviousLongitude, \
                                                                         VehicleLatitude, \
                                                                         VehicleLongitude)
        angle = self.calculate_position_angle(distance2, distance3, distance4)
      else:
        angle = -1

      ConcentricLayer1Item.Distance = distance2
      ConcentricLayer1Item.Angle = angle

    self.ConcentricLayer1PositionList.sort(key = lambda s: s.Distance)

    #print("[PONTEST][speedcamerad.py][recalculate_concentric_layer1()] len(self.ConcentricLayer1PositionList)=", len(self.ConcentricLayer1PositionList))
    #for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
    #  print("ConcentricLayer1Item=", ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance, ConcentricLayer1Item.Angle)

  def add_simulate_vihecle_data(self):
    print("[PONTEST][speedcamerad.py][add_simulate_vihecle_data()]")
    #Add vehicle tracke data
    VehiclePath = '/data/openpilot/selfdrive/vag/speedcamera_csv/GpsTrack.csv'
    try:
      fd = open(VehiclePath, 'r')
      rows = csv.reader(fd, delimiter=',')
      #print(time.ctime())
      for row in rows:
        self.VehicleMapPositionList.append(MapPosition(row[0], row[1], SpeedDirect.u))
      fd.close()
      print(time.ctime())
    except:
      print('ERROR: can not found ' + VehiclePath)
      exit(1)

    #for item in self.SpeedCameraMapPositionList:
    #  print("item=", item.Latitude, \
    #                 item.Longitude, \
    #                 item.Direct, \
    #                 item.SpeedLimit, \
    #                 item.RoadType)

    #Divide concentric layer
    #print("VehiclePosition=", self.VehicleMapPositionList[0].Latitude, self.VehicleMapPositionList[0].Longitude)
    #try:
    #  for MapPositionItem in self.SpeedCameraMapPositionList:
    #    #print("test", self.VehicleMapPositionList[0].Latitude, \
    #    #              self.VehicleMapPositionList[0].Longitude, \
    #    #              MapPositionItem.Latitude, \
    #    #              MapPositionItem.Longitude)
    #    distance1 = self.calculate_gps_radian_great_circle_distance(self.VehicleMapPositionList[0].Latitude, \
    #                                                                self.VehicleMapPositionList[0].Longitude, \
    #                                                                MapPositionItem.Latitude, \
    #                                                                MapPositionItem.Longitude)
    #    distance2 = self.calculate_gps_radian_haversine_formula_distance(self.VehicleMapPositionList[0].Latitude, \
    #                                                                     self.VehicleMapPositionList[0].Longitude, \
    #                                                                     MapPositionItem.Latitude, \
    #                                                                     MapPositionItem.Longitude)
    #    #print("SpeedCamPosition=", MapPositionItem.Latitude, MapPositionItem.Longitude, distance1, distance2)
    #except:
    #  print("Divide concentric layer fail")
    #  exit(1)

  def update_events(self, VehicleLatitude, VehicleLongitude, VehicleSpeed):
    #print("[PONTEST][speedcamerad.py][update_events()]", self.TestItemIndex, VehicleLatitude, VehicleLongitude, VehicleSpeed)

    if (len(self.ConcentricLayer1PositionList)==0):
      self.recalculate_concentric_layer_all(VehicleLatitude, VehicleLongitude, VehicleSpeed)
      self.recalculate_concentric_layer1(VehicleLatitude, VehicleLongitude, VehicleSpeed)
    else:
      if (self.TestItemIndex % 50) == 0:
        self.recalculate_concentric_layer_all(VehicleLatitude, VehicleLongitude, VehicleSpeed)
        self.recalculate_concentric_layer1(VehicleLatitude, VehicleLongitude, VehicleSpeed)
      else:
        self.recalculate_concentric_layer_all(VehicleLatitude, VehicleLongitude, VehicleSpeed)
        self.recalculate_concentric_layer1(VehicleLatitude, VehicleLongitude, VehicleSpeed)

    Item = 0
    for ConcentricLayer1Item in self.ConcentricLayer1PositionList:
      Item = Item + 1

      if self.SpeedCameraDetected < 3:
      #  print("[PONTEST][speedcamerad.py][update_events()] ConcentricLayer1Item=", self.TestItemIndex, self.VehicleMapPositionList[self.TestItemIndex].Latitude, self.VehicleMapPositionList[self.TestItemIndex].Longitude, ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance)
      #print("[PONTEST][speedcamerad.py][update_events()] len(self.oncentricLayer1PositionList)=", len(self.ConcentricLayer1PositionList))

        sc_send = messaging.new_message('speedCamera')
        sc_send.speedCamera.speedCameraMapPosition.latitude = float(ConcentricLayer1Item.Latitude)
        sc_send.speedCamera.speedCameraMapPosition.longitude = float(ConcentricLayer1Item.Longitude)
        sc_send.speedCamera.speedCameraMapPosition.direct = ConcentricLayer1Item.Direct
        sc_send.speedCamera.speedCameraMapPosition.speedLimitation = float(ConcentricLayer1Item.SpeedLimit)
        sc_send.speedCamera.speedCameraMapPosition.roadType = ConcentricLayer1Item.RoadType
        sc_send.speedCamera.speedCameraMapPosition.vehicleDistance = float(ConcentricLayer1Item.Distance)
        sc_send.speedCamera.speedCameraMapPosition.vehicleTrackAngle = float(ConcentricLayer1Item.Angle)
        sc_send.speedCamera.vehicleLatitude = float(VehicleLatitude)
        sc_send.speedCamera.vehicleLongitude = float(VehicleLongitude)
        sc_send.speedCamera.vehicleSpeed = float(VehicleSpeed)

        if (float(self.VehiclePreviousLatitude) < 180.0) and (float(self.VehiclePreviousLongitude) < 180.0):
          if(float(VehicleLatitude) > float(self.VehiclePreviousLatitude)):
            if(float(VehicleLongitude) > float(self.VehiclePreviousLongitude)):
              sc_send.speedCamera.vehicleDirect = SpeedDirect.ne
            else:
              sc_send.speedCamera.vehicleDirect = SpeedDirect.nw
          else:
            if(float(VehicleLongitude) > float(self.VehiclePreviousLongitude)):
              sc_send.speedCamera.vehicleDirect = SpeedDirect.se
            else:
              sc_send.speedCamera.vehicleDirect = SpeedDirect.sw
        else:
          sc_send.speedCamera.vehicleDirect = SpeedDirect.u

        if(float(ConcentricLayer1Item.Latitude) > float(VehicleLatitude)):
          if(float(ConcentricLayer1Item.Longitude) > float(VehicleLongitude)):
            self.VehicleToCameraDirect = SpeedDirect.ne
          else:
            self.VehicleToCameraDirect = SpeedDirect.nw
        else:
          if(float(ConcentricLayer1Item.Longitude) > float(VehicleLongitude)):
            self.VehicleToCameraDirect = SpeedDirect.se
          else:
            self.VehicleToCameraDirect = SpeedDirect.sw
        sc_send.speedCamera.vihicleToCameraDirect = self.VehicleToCameraDirect

        if (ConcentricLayer1Item.Direct == SpeedDirect.n and \
            (sc_send.speedCamera.vehicleDirect == SpeedDirect.ne or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.nw or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.n)) or \
           (ConcentricLayer1Item.Direct == SpeedDirect.s and \
            (sc_send.speedCamera.vehicleDirect == SpeedDirect.se or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.sw or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.s)) or \
           (ConcentricLayer1Item.Direct == SpeedDirect.e and \
            (sc_send.speedCamera.vehicleDirect == SpeedDirect.ne or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.se or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.e)) or \
           (ConcentricLayer1Item.Direct == SpeedDirect.w and \
            (sc_send.speedCamera.vehicleDirect == SpeedDirect.nw or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.sw or \
             sc_send.speedCamera.vehicleDirect == SpeedDirect.w)) or \
           ConcentricLayer1Item.Direct == SpeedDirect.d:
          VehicleDirectMatched = True
        else:
          VehicleDirectMatched = False

        #TODO: Need add 2nd and 3rd check
        if (ConcentricLayer1Item.Distance < self.VehiclePreviousSpeedCameraDistance) and \
           VehicleDirectMatched and \
           (((ConcentricLayer1Item.Distance < 0.5) and \
             (ConcentricLayer1Item.Distance > 0.1) and \
             (ConcentricLayer1Item.Angle < 10)) or \
            ((ConcentricLayer1Item.Distance < 0.1) and \
             self.SpeedCameraDetected)):
          self.SpeedCameraDetected = True
          sc_send.speedCamera.speedCameraDetected = True
        else:
          self.SpeedCameraDetected = False
          sc_send.speedCamera.speedCameraDetected = False

        print("SpeedCamItem=", self.TestItemIndex, VehicleLatitude, VehicleLongitude, \
                sc_send.speedCamera.vehicleDirect, sc_send.speedCamera.speedCameraDetected, \
                ConcentricLayer1Item.Latitude, ConcentricLayer1Item.Longitude, \
                ConcentricLayer1Item.Direct, ConcentricLayer1Item.SpeedLimit, \
                ConcentricLayer1Item.RoadType, ConcentricLayer1Item.Distance, \
                ConcentricLayer1Item.Angle, Item)
        self.pm.send('speedCamera', sc_send)

        self.VehiclePreviousLatitude = VehicleLatitude
        self.VehiclePreviousLongitude = VehicleLongitude
        self.VehiclePreviousSpeedCameraDistance = ConcentricLayer1Item.Distance
      else:
        break

  def speedcamerad_thread(self):
    print("[PONTEST][speedcamerad.py][speedcamerad_thread()]")

    #Simulate
    print("[PONTEST][speedcamerad.py][speedcamerad_thread()] Simulate")
    self.add_simulate_vihecle_data()
    while True:
      if self.TestItemIndex < 1900:
        #print("[PONTEST][speedcamerad.py][speedcamerad_thread()] self.TestItemIndex=", self.TestItemIndex)
        self.update_events(self.VehicleMapPositionList[self.TestItemIndex].Latitude, self.VehicleMapPositionList[self.TestItemIndex].Longitude, 50)
        self.TestItemIndex = self.TestItemIndex + 1
      else:
        print("[PONTEST][speedcamerad.py][speedcamerad_thread()] Simulate Reset")
        break;
      time.sleep(0.1) #100 ms

    #Real GPS data
    print("[PONTEST][speedcamerad.py][speedcamerad_thread()] Real GPS data")
    location_sock = messaging.sub_sock('gpsLocationExternal')
    while True:
      location = messaging.recv_sock(location_sock)
      if location:
        print("[PONTEST][speedcamerad.py][speedcamerad_thread()] location=", location.gpsLocationExternal.latitude, location.gpsLocationExternal.longitude, location.gpsLocationExternal.altitude, location.gpsLocationExternal.speed, location.gpsLocationExternal.timestamp, location.gpsLocationExternal.source)
        self.update_events(location.gpsLocationExternal.latitude, location.gpsLocationExternal.longitude, location.gpsLocationExternal.speed)

      time.sleep(0.1) #100 ms

    print("[PONTEST][speedcamerad.py][speedcamerad_thread()] end")


def main():
  print("[PONTEST][speedcamerad.py]Line:", sys._getframe().f_lineno)
  print("[PONTEST][speedcamerad.py][main()]")
  time.sleep(5)
  speedcamera = SpeedCamera()
  speedcamera.speedcamerad_thread()
  print("[PONTEST][speedcamerad.py][main()]")


if __name__ == "__main__":
  main()
