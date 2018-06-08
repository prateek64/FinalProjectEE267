#include "PoseTracker.h"
#include <Wire.h>

PoseTracker::PoseTracker(double alphaImuFilterIn, int baseStationModeIn, bool simulateLighthouseIn) :

  OrientationTracker(alphaImuFilterIn, false),
  lighthouse(),
  simulateLighthouse(simulateLighthouseIn),
  simulateLighthouseCounter(0),
  position{0, 0, -500},
  baseStationPitch(0),
  baseStationRoll(0),
  baseStationMode(baseStationModeIn),
  position2D{0,0,0,0, 0,0,0,0},
  clockTicks{0,0,0,0, 0,0,0,0},
  numPulseDetections{0,0,0,0, 0,0,0,0},
  pulseWidth{0,0,0,0,0,0,0,0}

  {

}

int PoseTracker::processLighthouse(float realXP, float realYP , float realF) {

  if (simulateLighthouse) {
  //if in simulation mode, get data from external file
    for (int i = 0; i < 8; i++) {
      clockTicks[i] = clockTicksData[(simulateLighthouseCounter*8 + i) % nLighthouseSamples];
      numPulseDetections[i] = 0;
    }

    //base station pitch/roll values remain the same throughout the simulation
    if (simulateLighthouseCounter == 0) {
      baseStationPitch = baseStationPitchSim;
      baseStationRoll = baseStationRollSim;
    }

    //data wraps around after end of array is reached
    simulateLighthouseCounter = (simulateLighthouseCounter + 1) % nLighthouseSamples;

    //slight delay to simulate delay between sensor readings (not exactly 120 Hz)
    delay(1);

  } else {
    //check data is available
    if (!lighthouse.readTimings(baseStationMode, clockTicks, numPulseDetections, pulseWidth,
      baseStationPitch, baseStationRoll)) {
      return -2;
    }

    //check that all diodes have detections
    for (int i = 0; i < 8; i++) {
      if (numPulseDetections[i] == 0) {
        return -1;
      }
    }
  }

  return updatePose(realXP, realYP , realF);

}

/**
 * TODO: see header file for documentation
 */
int PoseTracker::updatePose(float realXP, float realYP , float realF) {

//  float realXP = -0.399561;
//  float realYP = 1.3576;
//  float realF  = 2.0;

  
//  float realXP = 0.532554;
//  float realYP = 1.11969;
//  float realF  = 0.06;

//
//  float realXP = -0.116284;
//  float realYP = 0.729065;
//  float realF  =2.0;


 
  convertTicksTo2DPositions(clockTicks, position2D);



  for(int i = 0 ; i < 4 ; i++){

        position2D[2*i] = realF*position2D[2*i] + realXP;
        position2D[2*i + 1] = realF*position2D[2*i + 1] + realYP;

  }
    


  double A[8][8];
  formA(position2D, positionRef, A);

  double h[8];
  if (!solveForH(A, position2D, h)) {
    
    return 0;
  }


  double R[3][3];
  getRtFromH(h, R, position);

  quaternionHm = getQuaternionFromRotationMatrix(R);

  return 1;

}
