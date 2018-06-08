// Levenberg Marquardt Optimization 

#include "LevenbergMarquardt.h"
#include "MatrixMath.h"
#include "MatrixMath.cpp"
// #include "Utils.h"
#include <math.h> 


LevenbergMarquardtOpt::LevenbergMarquardtOpt(const float* photoDiodeLocations) {
	memcpy(photo2D, photoDiodeLocations, 8 * sizeof(float));
}

// LMparams Order ; [ThetaX , ThetaY , ThetaZ , TransX , TransY , TransZ , F , Xprincipal , Yprincipal]
void LevenbergMarquardtOpt::calcGOfX(const float *LMparams, float *h) {
	h[0] = LMparams[6]*(cos(LMparams[1])*cos(LMparams[2]) - sin(LMparams[0])*sin(LMparams[1])*sin(LMparams[2])) - LMparams[7]*cos(LMparams[0])*sin(LMparams[1]);
	h[1] = -LMparams[6]*cos(LMparams[0])*sin(LMparams[2]) + LMparams[7]*sin(LMparams[0]);
	h[2] = LMparams[6]*LMparams[3] + LMparams[7]*LMparams[5];
	h[3] = LMparams[6]*(cos(LMparams[1])*sin(LMparams[2]) + sin(LMparams[0])*sin(LMparams[1])*cos(LMparams[2])) - LMparams[8]*cos(LMparams[0])*sin(LMparams[1]);
	h[4] = LMparams[6]*cos(LMparams[0])*cos(LMparams[2]) + LMparams[8]*sin(LMparams[0]);
	h[5] = LMparams[6]*LMparams[4] + LMparams[8]*LMparams[5];
	h[6] = cos(LMparams[0])*sin(LMparams[1]);
	h[7] = -sin(LMparams[0]);
	h[8] = -LMparams[5];
}


// Output: f - 1x8 array of estimated x,y coordinates of all 4 photodiodes [x1,y1,x2,y2,x3,y3,x4,y4]
void LevenbergMarquardtOpt::calcFOfH(const float *h, float *f)  {

	for (int i = 0; i < 4; i++) {
		f[2*i] = (h[0]*photo2D[2*i] + h[1]*photo2D[2*i+1] + h[2])
						/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
		f[2*i+1] = (h[3]*photo2D[2*i] + h[4]*photo2D[2*i+1] + h[5])
						/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	}

}


// Output: JacobianG - 9x9 2D array jacobian of g(x)
void LevenbergMarquardtOpt::jacobianG(const float *LMparams, float JacobianG[][9])
{
	JacobianG[0][0] = -LMparams[6]*cos(LMparams[0])*sin(LMparams[1])*sin(LMparams[2]) + LMparams[7]*sin(LMparams[0])*sin(LMparams[1]);
	JacobianG[0][1] = -LMparams[6]*(sin(LMparams[1])*cos(LMparams[2]) - sin(LMparams[0])*cos(LMparams[1])*sin(LMparams[2])) - LMparams[7]*cos(LMparams[0])*cos(LMparams[1]);
	JacobianG[0][2] = -LMparams[6]*(cos(LMparams[1])*sin(LMparams[2]) + sin(LMparams[0])*sin(LMparams[1])*cos(LMparams[2]));
	JacobianG[0][3] = 0;
	JacobianG[0][4] = 0;
	JacobianG[0][5] = 0;
	JacobianG[0][6] = cos(LMparams[1])*cos(LMparams[2]) - sin(LMparams[0])*sin(LMparams[1])*sin(LMparams[2]);
	JacobianG[0][7] = cos(LMparams[0])*sin(LMparams[1]);
	JacobianG[0][8] = 0;

	JacobianG[1][0] = LMparams[6]*sin(LMparams[0])*sin(LMparams[2]) + LMparams[7]*cos(LMparams[0]);
	JacobianG[1][1] = 0;
	JacobianG[1][2] = -LMparams[6]*cos(LMparams[0])*cos(LMparams[2]);
	JacobianG[1][3] = 0;
	JacobianG[1][4] = 0;
	JacobianG[1][5] = 0;
	JacobianG[1][6] = -cos(LMparams[0])*sin(LMparams[2]);
	JacobianG[1][7] = sin(LMparams[0]);
	JacobianG[1][8] = 0;

	JacobianG[2][0] = 0;
	JacobianG[2][1] = 0;
	JacobianG[2][2] = 0;
	JacobianG[2][3] = LMparams[6];
	JacobianG[2][4] = 0;
	JacobianG[2][5] = LMparams[7];
	JacobianG[2][6] = LMparams[3];
	JacobianG[2][7] = LMparams[5];
	JacobianG[2][8] = 0;

	JacobianG[3][0] = LMparams[6]*cos(LMparams[0])*sin(LMparams[1])*cos(LMparams[2]) + LMparams[8]*sin(LMparams[0])*sin(LMparams[1]);
	JacobianG[3][1] = LMparams[6]*(-sin(LMparams[1])*sin(LMparams[2]) + sin(LMparams[0])*cos(LMparams[1])*cos(LMparams[2])) - LMparams[8]*cos(LMparams[0])*cos(LMparams[1]);
	JacobianG[3][2] = LMparams[6]*(cos(LMparams[1])*cos(LMparams[2]) - sin(LMparams[0])*sin(LMparams[1])*sin(LMparams[2]));
	JacobianG[3][3] = 0;
	JacobianG[3][4] = 0;
	JacobianG[3][5] = 0;
	JacobianG[3][6] = cos(LMparams[1])*sin(LMparams[2]) + sin(LMparams[0])*sin(LMparams[1])*cos(LMparams[2]);
	JacobianG[3][7] = 0;
	JacobianG[3][8] = - cos(LMparams[0])*sin(LMparams[1]);

	JacobianG[4][0] = -LMparams[6]*sin(LMparams[0])*cos(LMparams[2]) + LMparams[8]*cos(LMparams[0]);
	JacobianG[4][1] = 0;
	JacobianG[4][2] = -LMparams[6]*cos(LMparams[0])*sin(LMparams[2]);
	JacobianG[4][3] = 0;
	JacobianG[4][4] = 0;
	JacobianG[4][5] = 0;
	JacobianG[4][6] = cos(LMparams[0])*cos(LMparams[2]);
	JacobianG[4][7] = 0;
	JacobianG[4][8] = sin(LMparams[0]);

	JacobianG[5][0] = 0;
	JacobianG[5][1] = 0;
	JacobianG[5][2] = 0;
	JacobianG[5][3] = 0;
	JacobianG[5][4] = LMparams[6];
	JacobianG[5][5] = LMparams[8];
	JacobianG[5][6] = LMparams[4];
	JacobianG[5][7] = 0;
	JacobianG[5][8] = LMparams[5];

	JacobianG[6][0] = -sin(LMparams[0])*sin(LMparams[1]);
	JacobianG[6][1] = cos(LMparams[0])*cos(LMparams[1]);
	JacobianG[6][2] = 0;
	JacobianG[6][3] = 0;
	JacobianG[6][4] = 0;
	JacobianG[6][5] = 0;
	JacobianG[6][6] = 0;
	JacobianG[6][7] = 0;
	JacobianG[6][8] = 0;

	JacobianG[7][0] = -cos(LMparams[0]);
	JacobianG[7][1] = 0;
	JacobianG[7][2] = 0;
	JacobianG[7][3] = 0;
	JacobianG[7][4] = 0;
	JacobianG[7][5] = 0;
	JacobianG[7][6] = 0;
	JacobianG[7][7] = 0;
	JacobianG[7][8] = 0;

	JacobianG[8][0] = 0;
	JacobianG[8][1] = 0;
	JacobianG[8][2] = 0;
	JacobianG[8][3] = 0;
	JacobianG[8][4] = 0;
	JacobianG[8][5] = -1;
	JacobianG[8][6] = 0;
	JacobianG[8][7] = 0;
	JacobianG[8][8] = 0;

}


float LevenbergMarquardtOpt::sq(float x){
	return x*x;
}


void LevenbergMarquardtOpt::jacobianF(const float *h, float JacobianF[][9])
{
	for (int i = 0; i < 4; i++) {
	JacobianF[2*i][0] = photo2D[2*i]/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	JacobianF[2*i][1] = photo2D[2*i+1]/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	JacobianF[2*i][2] = 1/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	JacobianF[2*i][3] = 0;
	JacobianF[2*i][4] = 0;
	JacobianF[2*i][5] = 0;
	JacobianF[2*i][6] = -((h[0]*photo2D[2*i] + h[1]*photo2D[2*i+1] + h[2])/sq(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]))*photo2D[2*i];
	JacobianF[2*i][7] = -((h[0]*photo2D[2*i] + h[1]*photo2D[2*i+1] + h[2])/sq(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]))*photo2D[2*i+1];
	JacobianF[2*i][8] = -((h[0]*photo2D[2*i] + h[1]*photo2D[2*i+1] + h[2]))/sq(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);

	JacobianF[2*i+1][0] = 0;
	JacobianF[2*i+1][1] = 0;
	JacobianF[2*i+1][2] = 0;
	JacobianF[2*i+1][3] = photo2D[2*i]/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	JacobianF[2*i+1][4] = photo2D[2*i+1]/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	JacobianF[2*i+1][5] = 1/(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]);
	JacobianF[2*i+1][6] = -((h[3]*photo2D[2*i] + h[4]*photo2D[2*i+1] + h[5])/sq(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]))*photo2D[2*i];
	JacobianF[2*i+1][7] = -((h[3]*photo2D[2*i] + h[4]*photo2D[2*i+1] + h[5])/sq(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]))*photo2D[2*i+1];;
	JacobianF[2*i+1][8] = -((h[3]*photo2D[2*i] + h[4]*photo2D[2*i+1] + h[5])/sq(h[6]*photo2D[2*i] + h[7]*photo2D[2*i+1] + h[8]));
	}
}



// outputs: deltaX - 1x9 array of the step in which to move x
//					 [deltaThetaX , deltaThetaY , deltaThetaZ , deltaTransX , deltaTransY , deltaTransZ , deltaF , deltaXprincipal , deltaYprincipal ]
bool LevenbergMarquardtOpt::calcLMStep(const float* J, const float* f, const float* projection2D, float* deltaX) {
	int inv = 0;

	float JacobianAtT[9][8] = {};
	Matrix.Transpose((float*)J, 8, 9, (float*)JacobianAtT);

	float JacobianAtTJ[9][9] = {};
	Matrix.Multiply((float*)JacobianAtT, (float*)J, 9, 8, 9, (float*)JacobianAtTJ);


	float diag[9][9] = {};
	Matrix.MakeMatrixFromDiagonal((float*)JacobianAtTJ, 9, (float*)diag);

  	Matrix.Scale((float*)diag, 9, 9, lambda);

	float matToInvert[9][9] = {};
	Matrix.Add((float*)JacobianAtTJ, (float*)diag, 9, 9, (float*)matToInvert);

	if (! Matrix.Invert((float*)matToInvert, 9) ) {
		return false;
	}

	float leastSquaresAdagger[9][8] = {};
	Matrix.Multiply((float*)matToInvert, (float*)JacobianAtT, 9, 9, 8, (float*)leastSquaresAdagger);

	float leastSquaresb[8] = {};
	Matrix.Subtract((float*)projection2D, (float*)f, 8, 1, (float*)leastSquaresb);

	Matrix.Multiply((float*)leastSquaresAdagger, (float*)leastSquaresb, 9, 8, 1, deltaX);


	return true;


}


// Output : 1X9 array of // LMparams : [ThetaX , ThetaY , ThetaZ , TransX , TransY , TransZ , F , Xprincipal , Yprincipal]
bool LevenbergMarquardtOpt::estimateExtrinsicAndIntrinsicParameters(const float* projection2D,
        const float* poseEstimate) {
  // String label = "";
  bool successLM = true;

	for (int i = 0; i < maxNoIters; i++) {

		float h[9] = {};

		calcGOfX(poseEstimate, (float*)h);
		float f[8] = {};

		calcFOfH((const float*)h, (float*)f);
		float JacobianG[9][9] = {};

		jacobianG(poseEstimate, JacobianG);
		float JacobianF[8][9] = {};

		jacobianF((const float*)h, JacobianF);
		float J[8][9] = {};


		Matrix.Multiply((float*)JacobianF, (float*)JacobianG, 8, 9, 9, (float*)J);

		float deltaX[9] = {};

		if (! calcLMStep((const float*)J, (const float*)f, projection2D, (float*)deltaX)) {

			successLM = false;
			// Serial.println("not invertible!");
		}
		else {

			Matrix.Add((float*)poseEstimate, (float*)deltaX, 9, 1, (float*)poseEstimate);

		}


	}


	memcpy(estimatedExtrinsicAndIntrinsicParameters, poseEstimate, 9*sizeof(float));

	return successLM;
}



float LevenbergMarquardtOpt::calcResidual(const float *x, const float *projection2D)
{
	//float h[12];
	float h[9] = {};
	calcGOfX(x, h);
	float f[8] = {};
	calcFOfH(h, f);

	return (projection2D[0] - f[0]) * (projection2D[0] - f[0]) + (projection2D[1] - f[1]) * (projection2D[1] - f[1]) + (projection2D[2] - f[2]) * (projection2D[2] - f[2]) + (projection2D[3] - f[3]) * (projection2D[3] - f[3]) + (projection2D[4] - f[4]) * (projection2D[4] - f[4]) + (projection2D[5] - f[5]) * (projection2D[5] - f[5]) + (projection2D[6] - f[6]) * (projection2D[6] - f[6]) + (projection2D[7] - f[7]) * (projection2D[7] - f[7]);
}


