
/*

File that reads the photoDiode location Data from a text file and runs the Levenberg Marquardt Optimization technique to 
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

using namespace std;


int main(){

float photodiodeLocations[8] =	{-42.0, 25.0, 42.0, 25.0, 42.0, -25.0, -42.0, -25.0};

float projection2D[8];

// Array that holds the guess values from LM Optimization
float lmPoseGuess[9];

float A[8][8];

float hOut[8];

float Rout[3][3];

float initPos[8];
float Pos[8];

//Initial Intrinsic Parameters

float initF = 1.0;
float initXP = 0.3 , initYP = 0.4;

//Initial Rotation
float initRot[3];
//Initial Translation
float initTrans[3];

float residual = 0;
LevenbergMarquardtOpt lighthouseCalibrate(photodiodeLocations);

// Data File containing photodiode data 

// ifstream file("photoDiodeData.txt");
ifstream file("photoDiodeData_4thrun.txt");
for(int row = 0; row < 1; ++row) {
    std::string line;
    std::getline(file, line);
    if ( !file.good() ) 
        break;

    std::stringstream iss(line);

    for (int col = 0; col < 8; ++col) {
        std::string val;
        std::getline(iss, val, ',');
        if ( !iss.good() ) 
            break;

        std::stringstream convertor(val);
        convertor >> projection2D[col];

    }



}


for(int i = 0 ; i < 4 ; i++){


	initPos[2*i] = (projection2D[2*i] + initXP)/initF;
	initPos[2*i + 1] = (projection2D[2*i + 1] + initYP)/initF;

}


formA(initPos,photodiodeLocations,A);

solveForH(A , initPos , hOut);

getRtFromH(hOut , Rout , initTrans);


initRot[0] = 180*asin(Rout[2][1])/M_PI;
initRot[1] = 180*atan2(-Rout[2][0] , Rout[2][2])/M_PI;
initRot[2] = 180*atan2(-Rout[0][1] , Rout[1][1])/M_PI;


lmPoseGuess[0] = initRot[0];
lmPoseGuess[1] = initRot[1];
lmPoseGuess[2] = initRot[2];
lmPoseGuess[3] = initTrans[0];
lmPoseGuess[4] = initTrans[1];
lmPoseGuess[5] = initTrans[2];
lmPoseGuess[6] = initF;
lmPoseGuess[7] = initXP;
lmPoseGuess[8] = initYP;

// Change the for loop depending on number of data points you have in the file 

for(int row = 0; row < 6000; ++row) {
    std::string line;
    std::getline(file, line);
    if ( !file.good() ) 
        break;

    std::stringstream iss(line);

    for (int col = 0; col < 8; ++col) {

        std::string val;
        std::getline(iss, val, ',');
        if ( !iss.good() ) 
            break;

        std::stringstream convertor(val);
        convertor >> projection2D[col];

 	
 	}


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
 	
 	for(int i = 0 ; i < 4 ; i++) {


	projection2D[2*i] = (lmPoseGuess[6]*projection2D[2*i] + lmPoseGuess[7]);
	projection2D[2*i + 1] = (lmPoseGuess[6]*projection2D[2*i + 1] + lmPoseGuess[8]);

}
	
	// lmPoseGuess[6] = 1.0; 
 //    lmPoseGuess[7] = 0;
    // lmPoseGuess[8] = 0;
    lighthouseCalibrate.estimateExtrinsicAndIntrinsicParameters(projection2D, lmPoseGuess);
    // lmPoseGuess[6] = 1.0;
    // lmPoseGuess[7] = 0;
    // lmPoseGuess[8] = 0;

    residual += lighthouseCalibrate.calcResidual(lmPoseGuess ,projection2D);

    cout<<"Residual : "<<lighthouseCalibrate.calcResidual(lmPoseGuess ,projection2D)<<endl;
    cout<<"#################################"<<endl;
 	cout<<endl;


} 

cout<<"Average Residual : "<< residual/6000;


return 0;

}