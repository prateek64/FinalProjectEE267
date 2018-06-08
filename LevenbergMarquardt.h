// // Levenberg Marquardt Optimization Class

#ifndef _LevenbergMarquardt_h_
#define _LevenbergMarquardt_h_

class LevenbergMarquardtOpt
{
public:

	// Maximum Number of Iterations to run your LM optimization on
	static const int maxNoIters = 1000;

	//Damping Factor
	
	float residual[maxNoIters+1] = {};
	float estimatedExtrinsicAndIntrinsicParameters[9] = {};

	float photo2D[8] = {};

	float lambda = 0.1;

	// Constructor simply copies photoDiodeLocations into photo2D
	LevenbergMarquardtOpt(const float* photoDiodeLocations);
	~LevenbergMarquardtOpt(){};

	// Evaluates the function h=g(x) at x
	void calcGOfX(const float *LMparams, float *h);

	// Function that updates the parameters by running the LevenbergMarquardt Optimization 
	//Input - 2D photodiode location estimates and the current Parameter estimates
	bool estimateExtrinsicAndIntrinsicParameters(const float* projection2D, const float* positionEstimate);

	// Function to calculate the square of a function 
	float sq(float x);

	// Function that runs the update equation of the LM Optimization
	bool calcLMStep(const float* J, const float* f, const float* projection2D, float* deltaX);

	// Function to calculate the Jacobian of g(x) at the position x.
	void jacobianG(const float *LMparams, float JacobianG[][9]);

	// Function to calculate the Jacobian of h(x) at the position x.
	void jacobianF(const float *h, float JacobianF[][9]);

	// Function to calculate the function f=g(h) at h 
	void calcFOfH(const float *h, float *f);

	// Function to calculate residual
	float calcResidual(const float *LMparams, const float *pos2D);

};

#endif
