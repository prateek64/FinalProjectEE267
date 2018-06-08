/*

File that reads the simulated photoDiode location data and runs the Levenberg Marquardt Optimization technique to 
estimate the Intrinsic Parameters of the lighthouse

First an initial guess for the Focal length  and the principal points are given which run the homography method to get an 
initial guess for the pose . 

The initial guess array is then employed to run the Levenberg Marquardt optimization to get the optimal parameters

*/



#include "LevenbergMarquardt.h"
#include "LevenbergMarquardt.cpp"
#include "PoseMath.cpp"
#include "PoseMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include "simulatedLighthouseData.h"
using namespace std;


int main(){

float photodiodeLocations[8] =	{-42.0, 25.0, 42.0, 25.0, 42.0, -25.0, -42.0, -25.0};

float projection2D[8];

float lmPoseGuess[9];

float A[8][8];

float hOut[8];

float Rout[3][3];

float initialPos[8];

float initF = 2.5;
float initXP = 0.5 , initYP = 0.8;

float initRot[3];
float initTrans[3];
float initialTranslation[3];


float realXP = 7.0;
float realYP = 2.7;
float realF  = 1.0; 

uint32_t clockTicks[8];

LevenbergMarquardtOpt lighthouseCalibrate(photodiodeLocations);


for(int i = 0 ; i < 8 ; i++){


    clockTicks[i] = clockTicksData[i];


}

convertTicksTo2DPositions(clockTicks, projection2D);

for(int i = 0 ; i < 4 ; i++) {


	initialPos[2*i] = (projection2D[2*i] + realXP)/realF;
	initialPos[2*i + 1] = (projection2D[2*i + 1] + realYP)/realF;


}


formA(initialPos,photodiodeLocations,A);

solveForH(A , initialPos , hOut);

getRtFromH(hOut , Rout , initTrans);

// solveFor3DPosition(hOut , initialTranslation);

initRot[0] = asin(Rout[2][1]);
initRot[1] = atan2(-Rout[2][0] , Rout[2][2]);
initRot[2] = atan2(-Rout[0][1] , Rout[1][1]);


lmPoseGuess[0] = initRot[0];
lmPoseGuess[1] = initRot[1];
lmPoseGuess[2] = initRot[2];
lmPoseGuess[3] = initTrans[0];
lmPoseGuess[4] = initTrans[1];
lmPoseGuess[5] = initTrans[2];
lmPoseGuess[6] = initF;
lmPoseGuess[7] = initXP;
lmPoseGuess[8] = initYP;

cout<<"ThX : "<<lmPoseGuess[0]<<endl;
cout<<"ThY : "<<lmPoseGuess[1]<<endl;
cout<<"ThZ : "<<lmPoseGuess[2]<<endl;
cout<<"TrX : "<<lmPoseGuess[3]<<endl;
cout<<"TrY : "<<lmPoseGuess[4]<<endl;
cout<<"TrZ : "<<lmPoseGuess[5]<<endl;
cout<<"F : "<<lmPoseGuess[6]<<endl;
cout<<"XP : "<<lmPoseGuess[7]<<endl;
cout<<"YP : "<<lmPoseGuess[8]<<endl;
cout<<endl;
cout<<"#################################"<<endl;
cout<<endl;

for(int k = 0 ; k < 1770 ; k++){

    for(int i = 0 ; i < 8 ; i++){


        clockTicks[i] = clockTicksData[8*k + i];


    }


    convertTicksTo2DPositions(clockTicks, projection2D);
    cout<<"ThX : "<<lmPoseGuess[0]<<endl;
    cout<<"ThY : "<<lmPoseGuess[1]<<endl;
    cout<<"ThZ : "<<lmPoseGuess[2]<<endl;
    cout<<"TrX : "<<lmPoseGuess[3]<<endl;
    cout<<"TrY : "<<lmPoseGuess[4]<<endl;
    cout<<"TrZ : "<<lmPoseGuess[5]<<endl;
    cout<<"F : "<<lmPoseGuess[6]<<endl;
    cout<<"XP : "<<lmPoseGuess[7]<<endl;
    cout<<"YP : "<<lmPoseGuess[8]<<endl;
    cout<<endl;


   
    for(int i = 0 ; i < 8 ; i++){

        projection2D[2*i] = (projection2D[2*i] - realXP)/realF;
        projection2D[2*i + 1] = (projection2D[2*i + 1] - realYP)/realF;

    }



    lmPoseGuess[6] = 1.0;
    lmPoseGuess[7] = 0;
    lmPoseGuess[8] = 0;
    lighthouseCalibrate.estimateExtrinsicAndIntrinsicParameters(projection2D, lmPoseGuess);
    lmPoseGuess[6] = 1.0;
    lmPoseGuess[7] = 0;
    lmPoseGuess[8] = 0;
    

    float residual = lighthouseCalibrate.calcResidual(lmPoseGuess ,projection2D);

    cout<<"Residual : "<<residual<<endl;


    cout<<"#################################"<<endl;
    cout<<endl;

}

return 0;

}