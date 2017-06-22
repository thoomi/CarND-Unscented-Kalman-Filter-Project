#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_z_ = 3;
  n_aug_ = 7;
  n_sig_ = 2 * n_x_ + 1;
  n_sig_aug = 2 * n_aug_ + 1;
  lambda_ = 3;

  // initialize sigma point matrix with zeros
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sig_aug);

  // initialize weights
  weights_ = VectorXd(n_sig_aug);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int index = 1; index < n_sig_aug; ++index) {
    weights_(index) = 1 / (2 * (lambda_ + n_aug_));
  }

  z_pred_ = VectorXd::Zero(n_z_);
  Zsig_ = MatrixXd::Zero(n_z_, n_sig_aug);
  S_radar_ = MatrixXd::Zero(n_z_, n_z_);


  R_radar_ = MatrixXd(n_z_, n_z_);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;

  Q_radar_ = MatrixXd(2, 2);
  Q_radar_ << std_a_ * std_a_, 0,
              0, std_yawdd_ * std_yawdd_;


  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

/**
   * Initialize
   *
   */
void UKF::Initialize(MeasurementPackage first_meas_package) {
  if (first_meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    double range = first_meas_package.raw_measurements_(0);
    double phi = first_meas_package.raw_measurements_(1);

    x_ << range * cos(phi),
          range * sin(phi),
          0,
          0,
          0;
  }
  else if (first_meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    x_ << first_meas_package.raw_measurements_(0),
          first_meas_package.raw_measurements_(1),
          0,
          0,
          0;
  }

  previous_timestamp_ = first_meas_package.timestamp_;

  is_initialized_ = true;
}

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
    Initialize(meas_package);
    return;
  }

  double delta_t = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }

  std::cout << "x = " << x_ << std::endl;
  std::cout << "P = " << P_ << std::endl;
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
  MatrixXd Xsig_aug = GenerateSigmaPoints();
  Xsig_pred_ = PredictSigmaPoints(Xsig_aug, delta_t);

  // Predict state mean
  x_.fill(0.0);
  for (int index = 0; index < n_sig_aug; ++index) {
    x_ = x_ + weights_(index) * Xsig_pred_.col(index);
  }

  // Predict state covariance matrix
  P_.fill(0.0);
  for (int index = 0; index < n_sig_aug; ++index) {
    VectorXd x_diff = Xsig_pred_.col(index) - x_;

    x_diff(3) = NormalizeAngle(x_diff(3));

    P_ = P_ + weights_(index) * x_diff * x_diff.transpose();
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
  VectorXd z = VectorXd::Zero(2);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);


  VectorXd z_pred = H_laser_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
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
  PredictRadarMeasurment();

  VectorXd z = VectorXd(n_z_);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z_);

  //calculate cross correlation matrix
  for (int index = 0; index < n_sig_aug; ++index) {
    VectorXd x_diff = Xsig_pred_.col(index) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));

    VectorXd z_diff = Zsig_.col(index) - z_pred_;
    z_diff(1) = NormalizeAngle(z_diff(1));

    Tc = Tc + weights_(index) * (x_diff * z_diff.transpose());
  }

  //calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_z_, n_z_);
  K = Tc * S_radar_.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred_;
  z_diff(1) = NormalizeAngle(z_diff(1));

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_radar_ * K.transpose();
}


MatrixXd UKF::GenerateSigmaPoints() {
  MatrixXd Xsig = MatrixXd(n_aug_, n_sig_aug);

  // Augmented state vector
  VectorXd x_aug = VectorXd::Zero(n_aug_);
  x_aug.head(n_x_) = x_;

  // Augmented covariance matrix
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_radar_;

  // Square root of P_aug
  MatrixXd A = P_aug.llt().matrixL();

  // First row is the current state vector
  Xsig.col(0) = x_aug;

  // Generate points for each following column
  double c1 = sqrt(lambda_ + n_aug_);
  MatrixXd c2 = c1 * A;

  // Iterate state vector dimension times and calculate the sigma points
  for (int index = 0; index < n_aug_; ++index) {
    Xsig.col(index + 1) = x_aug + c2.col(index);
    Xsig.col(index + n_aug_ + 1) = x_aug - c2.col(index);
  }

  return Xsig;
}

MatrixXd UKF::PredictSigmaPoints(MatrixXd& rXsig_aug, double delta_t) {
  MatrixXd Xsig_pred = MatrixXd::Zero(n_x_, n_sig_aug);

  double dt2 = delta_t * delta_t;

  // Iterate over each sigma point column
  for (int index = 0; index < n_sig_aug; index++) {
    VectorXd x_sig = rXsig_aug.col(index);

    double px = x_sig(0);
    double py = x_sig(1);
    double v = x_sig(2);
    double yaw = x_sig(3);
    double yawd = x_sig(4);
    double noise_a = x_sig(5);
    double noise_yawdd = x_sig(6);

    // Calculate process noise vector
    VectorXd u = VectorXd(n_x_);
    u(0) = 0.5 * dt2 * cos(yaw) * noise_a;
    u(1) = 0.5 * dt2 * sin(yaw) * noise_a;
    u(2) = delta_t * noise_a;
    u(3) = 0.5 * dt2 * noise_yawdd;
    u(4) = delta_t * noise_yawdd;

    // Calculate predicted vector
    VectorXd x_pred = VectorXd(n_x_);
    if (yawd < 0.0001 && yawd > -0.0001) {
      x_pred(0) = v * cos(yaw) * delta_t;
      x_pred(1) = v * sin(yaw) * delta_t;
      x_pred(2) = 0;
      x_pred(3) = 0;
      x_pred(4) = 0;
    }
    else {
      x_pred(0) = v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      x_pred(1) = v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
      x_pred(2) = 0;
      x_pred(3) = yawd * delta_t;
      x_pred(4) = 0;
    }

    // Calculate final predicted sigma points
    Xsig_pred.col(index) = x_sig.head(n_x_) + x_pred + u;
  }

  return Xsig_pred;
}

void UKF::PredictRadarMeasurment() {
  z_pred_.setZero();
  S_radar_.setZero();

  //transform sigma points into measurement space
  for (int index = 0; index < n_sig_aug; index++) {
    VectorXd x_sig = Xsig_pred_.col(index);

    double px = x_sig(0);
    double py = x_sig(1);
    double v = x_sig(2);
    double psi = x_sig(3);

    // Calculate process noise vector
    VectorXd x_pred = VectorXd(n_z_);
    x_pred(0) = sqrt(px * px + py * py);
    x_pred(1) = atan2(py, px);
    x_pred(2) = (px * cos(psi) * v + py * sin(psi) * v) / x_pred(0);

    // Calculate final predicted sigma points
    Zsig_.col(index) = x_pred;
  }

  //calculate mean predicted measurement
  for (int index = 0; index < n_sig_aug; ++index) {
    z_pred_ = z_pred_ + weights_(index) * Zsig_.col(index);
  }

  //predict state covariance matrix
  for (int index = 0; index < n_sig_aug; ++index) {
    VectorXd c1 = Zsig_.col(index) - z_pred_;
    MatrixXd c2 = c1 * c1.transpose();

    S_radar_ = S_radar_ + weights_(index) * c2;
  }

  S_radar_ = S_radar_ + R_radar_;
}

double UKF::NormalizeAngle(double angle) {
  double result = angle;

  while (result > M_PI) result -= 2.0 * M_PI;
  while (result < -M_PI) result += 2.0 * M_PI;

  return result;
}