#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.setIdentity();

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.8;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  /**
  TODO:
 Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;

  //State dimension
  n_x_ = 5;

  //Augmented state dimension
  n_aug_ = 7;

  //Define spreading parameter
  lambda_ = 3 - n_aug_;

  // Matrix to hold sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  //Vector for weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++)
  { //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // Noise matrices
  R_radar = MatrixXd(3,3);
  R_laser = MatrixXd(2,2);

  //use EKF to update lidar
  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;
  R_laser << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;
  R_radar << std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;

  // Start time
  time_us_ = 0;

  // Set NIS
  NIS_radar_ = 0;
  NIS_laser_ = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  //Initialization//
  // skip processing if the both sensors are ignored
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
      (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_))
  {
    if (!is_initialized_)
    {
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
      {
        double rho = meas_package.raw_measurements_[0];
        double phi = meas_package.raw_measurements_[1];
        double rho_dot = meas_package.raw_measurements_[2];

        //coordinates transfer and update x_
        x_(0) = rho * cos(phi);
        x_(1) = rho * sin(phi);

      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
      {
        x_(0) = meas_package.raw_measurements_[0];
        x_(1) = meas_package.raw_measurements_[1];
      }
      // done initializing, no need to predict or update
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
  }

  //Prediction
  // Calculate delta_t, store current time for future
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict
  Prediction(delta_t);

  //Update//
  // Measurement updates
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create Augmented Sigma Points
  // creae augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;

  //creae augmented convariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L_ = P_aug_.llt().matrixL();

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug_.fill(0.0);

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
  }

  //Predict sigma points

  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    //extract values for better readability
    const double p_x = Xsig_aug_(0,i);
    const double p_y = Xsig_aug_(1,i);
    const double v = Xsig_aug_(2,i);
    const double yaw = Xsig_aug_(3,i);
    const double yawd = Xsig_aug_(4,i);
    const double nu_a = Xsig_aug_(5,i);
    const double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    const double check_zero = 0.001;
    if (fabs(yawd) > check_zero) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //extract measurement as VectgorXd
  VectorXd z_ = meas_package.raw_measurements_;

  VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z_ - z_pred;
	MatrixXd Ht = H_laser_.transpose();
  MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_laser_ * PHt + R_laser;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;
  NIS_laser_ = y.transpose() * Si * y;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //extract measurement as VectorXd
  VectorXd z = meas_package.raw_measurements_;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig_.fill(0.0);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    if ((p_x == 0.0) && (p_y == 0.0)) {
      // handle atan2(0,0)
      Zsig_(1, i) = 0.0;
    } else{
      Zsig_(1, i) = atan2(p_y, p_x);
    }                                                //phi
    const double epsilon = 0.00001;
    Zsig_(2, i) = (p_x*v1 + p_y*v2) / std::max(epsilon, Zsig_(0, i));   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }


  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);

  //calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    if (z_diff(1) > M_PI) {
      z_diff(1) -= 2. * M_PI;
    } else if (z_diff(1) < - M_PI) {
      z_diff(1) += 2. * M_PI;
    }
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add R (noise) to S
  S = S + R_radar;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  //calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //normalize angles
    if (x_diff(3) > M_PI) {
      x_diff(3) -= 2. * M_PI;
    } else if (x_diff(3) < -M_PI) {
      x_diff(3) += 2. * M_PI;
    }
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    //normalize angles
    if (z_diff(1) > M_PI) {
      z_diff(1) -= 2. * M_PI;
    } else if (z_diff(1) < -M_PI) {
      z_diff(1) += 2. * M_PI;
    }
    Tc += weights_(i) * x_diff * z_diff.transpose();

  }

  // residual
  VectorXd z_diff = z - z_pred_;

  //normalize angles
  if (z_diff(1) > M_PI) {
    z_diff(1) -= 2. * M_PI;
  } else if (z_diff(1) < -M_PI) {
    z_diff(1) += 2. * M_PI;
  }

  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
}
