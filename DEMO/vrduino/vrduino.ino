/**
 *  This is the main file for pose tracking with the VRduino.
 */

#include <Wire.h>
#include "TestPose.h"
#include "PoseTracker.h"
#include "InputCapture.h"

//complementary filter value [0,1].
//1: ignore acc tilt, 0: use all acc tilt
double alphaImuFilter = 0.99;

float DL = 1.0;
//get simulated lighthouse timings (to test without physical lighthouse)
bool simulateLighthouse = false;

//if test is true, then run tests in Test.cpp and exit
bool test = false;

//mode of base station
//0:A, 1:B, 2: C
const int A = 0;
const int B = 1;
const int C = 2;
int baseStationMode = C;

float realXP = -0.004393;
float realYP = 0.426495;
float realF  =  0.979;

//if true, measure the imu bias on start
bool measureImuBias = true;

//if measureImuBias is false, set the imu bias to the following
double imuBias[3] = {0, 0, 0};

PoseTracker tracker(alphaImuFilter, baseStationMode, simulateLighthouse);

void setup() {

  Serial.begin(115200);
  if (test) {

    delay(500);
    testPoseMain();
    return;

  }

  tracker.initImu();

  if (measureImuBias) {

    tracker.measureImuBiasVariance();

  } else {

    tracker.setImuBias(imuBias);

  }

}

void loop() {

  if (test) {

    return;

  }

  if (Serial.available()) {

    int byteRead = Serial.read();
    int desiredMode  = byteRead - 48;

    if(byteRead == 'q') {

        realXP = -0.004393;
        realYP = 0.426495;
        realF  =  0.979;    
        DL = 0.41;

        
//       Serial.printf("PA ");
//       Serial.printf("%.3f\n" , ALp);
    }

    else if(byteRead == 'w') {

        realXP = -0.004393;
        realYP = 0.426495;
        realF  =  0.979;
        DL = 1.0;

    }

    if (desiredMode >= 0 && desiredMode <= 2) {

      tracker.setMode(desiredMode);

    } else if (byteRead == 'r') {

      //reset orientation tracking
      tracker.resetOrientation();

    } else if (byteRead == 'b') {

      //remeasure bias
      tracker.measureImuBiasVariance();

    }

  }


  //variables determining the success of imu and lighthouse tracking
  bool imuTrack = false;
  int hmTrack = -2;

  imuTrack = tracker.processImu();
  hmTrack = tracker.processLighthouse(realXP , realYP , realF);

  //get values from tracker
  double pitch = tracker.getBaseStationPitch();
  double roll = tracker.getBaseStationRoll();
  int mode = tracker.getBaseStationMode();
  const unsigned long * numPulseDetections = tracker.getNumPulseDetections();
  const double * position = tracker.getPosition();
  const double * position2D = tracker.getPosition2D();
  const Quaternion& quaternionComp = tracker.getQuaternionComp();
  const Quaternion& quaternionHm = tracker.getQuaternionHm();


  if (hmTrack > -2) {
    // base station data available

    //print base station data
    Serial.printf("BS ");
    Serial.printf("%.3f %.3f %d\n", pitch, roll, mode);

    // print num sweep pulse detections for each axis of each photodiode
    // order is sensor0x, sensor0y, ... sensor3x, sensor3y
    // should normally be 1 1 1 1 1 1 1 1
    // could be more than 2 if there are inter-reflections,
    // or 0 if there are occlusions
    Serial.printf("NP ");
    for (int i = 0; i < 8; i++) {
      Serial.printf("%lu ", numPulseDetections[i]);
    }
    Serial.println();

  }

  if (hmTrack == 1 ) {

    //print xyz position
    Serial.printf("PS %.3f %.3f %.3f %.3f\n",
      position[0], position[1], position[2] , DL);

    //print quaternion from homography
    Serial.printf("QH %.3f %.3f %.3f %.3f\n",
      quaternionHm.q[0], quaternionHm.q[1],
      quaternionHm.q[2], quaternionHm.q[3]);

    //print 2D positions of photodiodes for visualization
    Serial.printf("PD ");
    for (int i = 0; i < 8; i++) {
      Serial.printf("%.3f ", position2D[i]);
    }
    Serial.println();

  }

  if (imuTrack == 1) {

  //print quaternion from imu
    Serial.printf("QC %.3f %.3f %.3f %.3f\n",
      quaternionComp.q[0], quaternionComp.q[1],
      quaternionComp.q[2], quaternionComp.q[3]);

  }

}
