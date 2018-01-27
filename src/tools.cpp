#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    if(estimations.size() == ground_truth.size()) {
        //accumulate squared residuals
        for(int i=0; i < estimations.size(); ++i){
            // ... your code here
            VectorXd r = estimations[i] - ground_truth[i];
            r = r.array()*r.array();
            rmse += r;
        }
        
        //calculate the mean
        rmse = rmse/estimations.size();
        //calculate the squared root
        rmse = rmse.array().sqrt();
        return rmse;
    }
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    MatrixXd Hj(3,4);
    //recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    //check division by zero
    if(fabs(px+py)<0.0001){
        cout << "Can not divid by Zero";
        return Hj;
    };
    
    //compute the Jacobian matrix
    Hj << px/sqrt(pow(px,2.0)+pow(py,2.0)), py/sqrt(pow(px,2.0)+pow(py,2.0)), 0, 0,
    -py/(pow(px,2.0)+pow(py,2.0)), px/(pow(px,2.0)+pow(py,2.0)), 0, 0,
    py*(vx*py-vy*px)/pow(pow(px,2.0)+pow(py,2.0),(3/2)), px*(vy*px-vx*py)/pow(pow(px,2.0)+pow(py,2.0),(3/2)),
    px/sqrt(pow(px,2.0)+pow(py,2.0)), py/sqrt(pow(px,2.0)+pow(py,2.0));
    
    return Hj;
}
