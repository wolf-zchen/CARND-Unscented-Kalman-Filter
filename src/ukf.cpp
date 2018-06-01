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

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.4;
  
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
  //State dimension
  n_x_ = 5;

  //Augmented state dimension
  n_aug_ = 7;

  //Define spreading parameter
  lambda_ = 0;
 
  // Matrix to hold sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //Vector for weights
  weights_ = VectorXd(2*n_aug_+1);

  // Noise matrices
  R_radar = MatrixXd(3,3);
  R_laser = MatrixXd(2,2);
  
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
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      //coordinates transfer and update x_
      // Middle value for 'v' can be tuned
      x_ << rho * cos(phi), rho * sin(phi), 3, rho_dot * cos(phi), rho_dot * sin(phi);

      //state covariance matrix
      P_ << std_radr_*std_radr_, 0, 0, 0, 0,
            0, std_radr_*std_radr_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, std_radphi_, 0,
            0, 0, 0, 0, std_radphi_;
      
      // Create R for update noise later
      R_radar << std_radr_*std_radr_, 0, 0,
                 0, std_radphi_*std_radphi_, 0,
                 0, 0, std_radrd_*std_radrd_;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 3,0, 0;

      //state covariance matrix
      P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
            0, std_laspy_*std_laspy_, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;
      
      // Create R for update noise later
      R_laser << std_laspx_*std_laspx_, 0,
                 0, std_laspy_*std_laspy_;

    }
    // done initializing, no need to predict or update
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // Calculate delta_t, store current time for future
  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict
  Prediction(delta_t);

  // Measurement updates
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
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

  //Sigma point prediction
  // Define spreading parameter
  lambda_ = 3 - n_x_;

  //set sigma point matrix
  MatrixXd Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A_ = P_.llt().matrixL();
  
  //calculate sigma points
  Xsig_.col(0)  = x_;
  for (int i = 0; i < n_x_; i++)
  {
    Xsig_.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A_.col(i);
    Xsig_.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A_.col(i);
  }

  // creae augmented mean vector
  VectorXd x_aug_ = VectorXd(n_aug_);

  //create augmented tate convariance
  MatrixXd P_aug_ = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create spreading parameter
  lambda_ = 3 - n_aug_;
  //create augmented mean state
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;
  
  //creae augmented convariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L_ = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
  }

  //predict sigma points
  //create vector for predicted state
  VectorXd x_pred_ = VectorXd(n_x_);

  //create convariance matrix for prediction
  MatrixXd P_pred_ = MatrixXd(n_x_, n_x_);

  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
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

  //set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) 
  { //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  //predicted state mean
  x_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //iterate over sigma points
    x_pred_ = x_pred_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_pred_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    P_pred_ = P_pred_ + weights_(i) * x_diff_ * x_diff_.transpose() ;
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
  //set measurement dimension, lidar measur4e px and py
  int n_z_ = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2*n_aug_+1);

  //create mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);

  //create predicted measurement convariance matrix
  MatrixXd S = MatrixXd(n_z_, n_z_);
  Zsig_.fill(0.0);
  z_pred_.fill(0.0);
  S.fill(0.0);

  for (int i=0; i<2*n_aug_+1; i++)
  {
  VectorXd state = Xsig_pred_.col(i);
  double px = state(0);
  double py = state(1);
  Zsig_.col(i) << px, py;

  //calculate mean predicted measurement
  z_pred_ += weights_(i) * Zsig_.col(i);
  }

  //calculate measurement covariance matrix S
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  // Add R (noise) to S
  S += R_laser;
  
  //create vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_);
  
  double meas_px = meas_package.raw_measurements_(0);
  double meas_py = meas_package.raw_measurements_(1);
  
  z << meas_px,
       meas_py;
  
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

    Tc += weights_(i) * x_diff * z_diff.transpose();

  }
  
  // residual
  VectorXd z_diff = z - z_pred_;
  
  //calculate NIS
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //update state mean and covariance matrix
  x_ += K*z_diff;
  P_ -= K*S*K.transpose();
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
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  
  //mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  
  Zsig_.fill(0.0);
  z_pred_.fill(0.0);
  S.fill(0.0);
  double rho = 0;
  double phi = 0;
  double rho_d = 0;
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    //transform sigma points into measurement space
    VectorXd state = Xsig_pred_.col(i);
    double px = state(0);
    double py = state(1);
    double v = state(2);
    double yaw = state(3);
    double yaw_d = state(4);
    
    rho = sqrt(px*px+py*py);
    phi = atan2(py,px);
    rho_d = (px*cos(yaw)*v+py*sin(yaw)*v) / rho;
    
    Zsig_.col(i) << rho, phi, rho_d;
    
    //calculate mean predicted measurement
    z_pred_ += weights_(i) * Zsig_.col(i);
  }
  
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
  S += R_radar;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_);
  
  double meas_rho = meas_package.raw_measurements_(0);
  double meas_phi = meas_package.raw_measurements_(1);
  double meas_rhod = meas_package.raw_measurements_(2);
  
  z << meas_rho, meas_phi, meas_rhod;
  
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
