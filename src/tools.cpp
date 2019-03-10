#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
   VectorXd rmse(4);
   rmse << 0,0,0,0;
   if (estimations.size() == 0 || estimations.size() != ground_truth.size())
   {
       std::cout << "invalid size of estimations and ground_truth vector" << std::endl;
       return rmse;
   }

   for (int i = 0; i < estimations.size(); i++)
   {
       VectorXd residual = estimations[i] - ground_truth[i];
       residual = residual.array() * residual.array();
       rmse +=  residual;
   }    
   rmse = rmse / estimations.size();
   rmse = rmse.array().sqrt();
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
     MatrixXd Hj(3,4);
     
     // Set default value to zero;
    Hj << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);


    double var = (px * px) + (py * py);
    double var_sqrt = sqrt(var);
    double Var_rt_3_2 = (var * var_sqrt);

    // check division by zero
    if (fabs(var) < 0.0001)
    {
        std::cout << "Divide by zero in CalculateJacobian" << std::endl;
        getchar();
        return Hj;
    }

    // compute the Jacobian matrix
    Hj << (px / var_sqrt), (py / var_sqrt), 0, 0,
      -(py / var), (px / var), 0, 0,
      (py * (vx * py - vy * px) / Var_rt_3_2), (px*(px*vy - py*vx)/Var_rt_3_2), px/var_sqrt, py/var_sqrt;

    return Hj;
}
