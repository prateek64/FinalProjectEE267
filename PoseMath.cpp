#include "PoseMath.h"

#include <math.h>
/**
 * TODO: see header file for documentation
 */


void convertTicksTo2DPositions(uint32_t clockTicks[8], float pos2D[8]) {
  //use variable CLOCKS_PER_SECOND defined in PoseMath.h
  //for number of clock ticks a second

	float ang[8];

	for(int i = 0 ; i < 8 ; i++){

		if(i%2 == 0){

			ang[i] = 90.0 - double(clockTicks[i])*60.0*360.0/48000000;

		}

		else {

			ang[i] = -90.0 + double(clockTicks[i])*60.0*360.0/48000000;
		}

		pos2D[i] = tan(2*M_PI*ang[i]/360.0);
    // Serial.write(2);

	}

    

}

/**
 * TODO: see header file for documentation
 */
void formA(float pos2D[8], float posRef[8], float Aout[8][8]) {

	for (int i = 0; i < 4; i++) {

        Aout[i*2][0] = posRef[i*2];
        Aout[i*2][1] = posRef[i*2+1];
        Aout[i*2][2] = 1.0;
        Aout[i*2][3] = 0.0;
        Aout[i*2][4] = 0.0;
        Aout[i*2][5] = 0.0;
        Aout[i*2][6] = -posRef[i*2]*pos2D[i*2];
        Aout[i*2][7] = -posRef[i*2+1]*pos2D[i*2];
        Aout[i*2+1][0] = 0.0;
        Aout[i*2+1][1] = 0.0;
        Aout[i*2+1][2] = 0.0;
        Aout[i*2+1][3] = posRef[i*2];
        Aout[i*2+1][4] = posRef[i*2+1];
        Aout[i*2+1][5] = 1.0;
        Aout[i*2+1][6] = -posRef[i*2]*pos2D[i*2+1];
        Aout[i*2+1][7] = -posRef[i*2+1]*pos2D[i*2+1];

    }


}


/**
 * TODO: see header file for documentation
 */
bool solveForH(float A[8][8], float b[8], float hOut[8]) {
  //use Matrix Math library for matrix operations
  //example:
  //int inv = Matrix.Invert((double*)A, 8);
  //if inverse fails (Invert returns 0), return false 
  
  int inv =  Matrix.Invert((float*)A, 8);
    if(inv == 0)
        return false;
    Matrix.Multiply((float*)A, (float*)b, 8, 8, 1, (float*) hOut);
    return true;


}


/**
 * TODO: see header file for documentation
 */

float sq(float x){
	return x*x;
}

void getRtFromH(float h[8], float ROut[3][3], float pos3DOut[3]) {

 double lengthR1 = sqrt( h[0]*h[0] + h[3]*h[3] + h[6]*h[6] );
  double lengthR2 = sqrt( h[1]*h[1] + h[4]*h[4] + h[7]*h[7] );

  // normalization factor for position
  double normalizationFactor = (lengthR1 + lengthR2) / 2.0;

  // Translation vector:
  pos3DOut[0] = h[2]/normalizationFactor;
  pos3DOut[1] = h[5]/normalizationFactor;
  pos3DOut[2] = -1.0/normalizationFactor;

  // get rotation
  // get first column of rotation matrix and normalize
  double r1[3] = { h[0]/lengthR1, h[3]/lengthR1, -h[6]/lengthR1 };


  // get second column
  double dotProd = r1[0]*h[1] + r1[1]*h[4] - r1[2]*h[7];

  double r2[3] = {
    h[1] - dotProd*r1[0],
    h[4] - dotProd*r1[1],
    -h[7] - dotProd*r1[2] };

  double normr2 = sqrt( r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2] );
  r2[0] = r2[0] / normr2;
  r2[1] = r2[1] / normr2;
  r2[2] = r2[2] / normr2;

  // get 3rd column that is orthogonal to both of the others
  double r3[3] = {
    r1[1]*r2[2] - r1[2]*r2[1],
    r1[2]*r2[0] - r1[0]*r2[2],
    r1[0]*r2[1] - r1[1]*r2[0] };

  // populate rotation matrix in row-major order
  for (int i = 0; i < 3; i++) {
    ROut[i][0] = r1[i];
    ROut[i][1] = r2[i];
    ROut[i][2] = r3[i];
  }
}


float solveFor3DPosition(float h[8] , float position3D[3]) {
  float normalizationFactor = 2.0/( sqrt( sq(h[0]) + sq(h[3]) + sq(h[6]) ) +  sqrt( sq(h[1]) + sq(h[4]) + sq(h[7]))  );

  float origin_of_board[3] = {0.0, 0.0, 1.0};
  float h_extended[3][3] = {h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], 1.0};
  int m = sizeof(h_extended)/sizeof(h_extended[0]); //dimension of h

  Matrix.Scale((float*)h_extended, m, m, normalizationFactor);
  Matrix.Multiply((float*)h_extended, (float*)origin_of_board, m, m, 1, (float*)position3D);
  position3D[2] = -1.0*position3D[2];
  return normalizationFactor;
}

/**
 * TODO: see header file for documentation
 */
// Quaternion getQuaternionFromRotationMatrix(double R[3][3]) {


//   double q0 = sqrt(1.0 + R[0][0] + R[1][1] + R[2][2])/2.0;

//   double qx = (R[2][1] - R[1][2])/(4.0*q0);
//   double qy = (R[0][2] - R[2][0])/(4.0*q0);
//   double qz = (R[1][0] - R[0][1])/(4.0*q0);

//   return Quaternion(q0 , qx , qy , qz);

// }
