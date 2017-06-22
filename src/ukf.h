#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  double previous_timestamp_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* Predicted radar measurement state
  VectorXd z_pred_;

  ///* Predicted radar measurement sigma points
  MatrixXd Zsig_;

  ///* Radar measurment noise matrix
  MatrixXd R_radar_;

  ///* Radar process noise matrix
  MatrixXd Q_radar_;

  ///* Predicted radar measurement covariance matrix
  MatrixXd S_radar_;

  ///* Laser process transformation matrix
  MatrixXd H_laser_;

  ///* Laser noise matrix
  MatrixXd R_laser_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Radar measurement state dimension
  int n_z_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point vector dimension without augmentation
  int n_sig_;

  ///* sigma point vector dimension with augmentation
  int n_sig_aug;

  ///* Sigma point spreading parameter
  double lambda_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Initialize
   *
   */
  void Initialize(MeasurementPackage first_meas_package);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

private:

  MatrixXd GenerateSigmaPoints();
  MatrixXd PredictSigmaPoints(MatrixXd& rXsig_aug, double delta_t);

  void PredictRadarMeasurment();

private:

  double NormalizeAngle(double angle);
};

#endif /* UKF_H */
