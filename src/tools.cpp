#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
	// re-used from EKF...
	
  	VectorXd rmse(4);
	rmse << 0,0,0,0;

    if(estimations.size() == 0 || estimations.size() != ground_truth.size())
    {
        cout << "error" << std::endl;
        return rmse;
    }

    VectorXd residual(4);
    residual << 0,0,0,0;
    
	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
        // ... your code here
        VectorXd sample = estimations[i] - ground_truth[i];
        //cout << "sample" << endl << sample << endl;
        sample = pow(sample.array(), 2);
        //cout << "sample2" << endl << sample << endl;
		residual += sample;
		//cout << "residual" << endl << residual << endl;
	}

	//calculate the mean
	//cout << "est size " << estimations.size() << endl;
	//cout << "est inv " << 1.0/estimations.size() << endl;
    residual = (1.0/estimations.size()) * residual.array();
    //cout << "residual mean" << endl << residual << endl;

	//calculate the squared root
	rmse = residual.array().sqrt();

	//return the result
	return rmse;
}