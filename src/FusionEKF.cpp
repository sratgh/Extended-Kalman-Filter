#include "FusionEKF.h"
#include <iostream>
#include <Eigen/Dense>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {

  // Initialize state vector
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1, 1, 1, 1;

  // Set is_initialized_ to false in the beginning
  is_initialized_ = false;

  // Set previous_timestamp_ to 0 in the beginning
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Initizalize measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Initizalize measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
   // Initialize Covariance P-Matrix
   ekf_.P_ = MatrixXd(4,4);
   ekf_.P_ << 1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;
   // initialize H_laser matrix to map the state vector x into the measurement space
   H_laser_ << 1, 0, 0, 0,
               0, 1, 0, 0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

/**
 * Process Measurement function
 * @param measurement_pack
 */
void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates
      //         and initialize state.
      // Extract raw measurements
      float roh = measurement_pack.raw_measurements_[0];      // range
      float phi = measurement_pack.raw_measurements_[1];      // bearing
      float roh_dot = measurement_pack.raw_measurements_[2];  // velocity of rho

      // convert from polar to cartesian coordinates
      float x = roh * cos(phi);
      float y = roh * sin(phi);
      float v_x = roh_dot * cos(phi);
      float v_y = roh_dot + sin(phi);

      // store values in state vector
      ekf_.x_ << x, y, v_x, v_y;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.

      // Initialize state with raw measurements directly, since already cartesian
      // no velocity in case of laser measurement
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
    }

    // Set previous timestamp for calculation of the time
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   */
  // Calculate the dt from the two timestamps
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Set the F-Matrix
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // Set the acceleration noise / process noise
  float noise_ax = 5.0;
  float noise_ay = 5.0;

  // Precalculate some of the values, needed for the completing the Q-Matrix
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  float dt_4_4 = dt_4 / 4;
  float dt_3_2 = dt_3 / 2;

  // Set the Q-Matrix with the calculated values
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
             0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
             dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
             0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;

  // Call the prediction function
  ekf_.Predict();

  /**
   * Update
   */

  /**
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
