#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
/*
 * Please note that the Eigen library does not initialize
 * VectorXd or MatrixXd objects with zeros upon creation.
 */

 /**
  * Constructor
  */
KalmanFilter::KalmanFilter() {

}

/**
 * Destructor
 */
KalmanFilter::~KalmanFilter() {}

/**
 * Predict function
 */
void KalmanFilter::Predict() {
  /**
   * predict the state
   */
   // Update state vector with model prediction
   // The model is a linear model, assuming a constant velocity
   x_ = F_ * x_;
   // Updating the covariance matrix
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

/**
 * Update Kalman filter function
 */
void KalmanFilter::Update(const VectorXd &z) {
  /**
   * update the state by using Kalman Filter equations
   */
   // Update with the measurement
   VectorXd y = z - H_ * x_;

   // Call the common update function (for radar and laser)
   update_(y);
}

/**
 * Update Extended Kalman filter function
 */
void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * update the state by using Extended Kalman Filter equations
   */
   float px = x_(0);
   float py = x_(1);
   float vx = x_(2);
   float vy = x_(3);

   // convert from cartesian to polar
   float rho = sqrt(px*px + py*py);
   float phi = atan2(py, px);
   float rho_dot = (px*vx + py*vy) / rho;
   VectorXd h_x = VectorXd(3);
   h_x << rho, phi, rho_dot;

   VectorXd y = z - h_x;
   // normalize the angle between -pi to pi
   while ( y(1) > M_PI || y(1) < -M_PI ) {
       if ( y(1) > M_PI ) {
         y(1) -= M_PI;
       } else {
         y(1) += M_PI;
       }
     }
   // Call the common update function (for radar and laser)
   update_(y);
}

/**
 * Common update function for both extended and normal kalman filter
 */
void KalmanFilter::update_(const VectorXd &y)
{
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
