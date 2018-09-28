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
  P_.fill(0.0);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.4;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.52;
  
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
  R_radar_ = MatrixXd(3,3);
  R_radar_ <<    std_radr_*std_radr_, 0, 0,
                 0, std_radphi_*std_radphi_, 0,
                 0, 0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  time_us_ = 0;
  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  is_initialized_ = false;
  P_ << .3, 0, 0,  0,  0,
        0, .5, 0,  0,  0,
        0,  0, 4,  0,  0,
        0,  0, 0, .1,  0,
        0,  0, 0,  0, .1;

  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  
  for (int i=1; i<2*n_aug_+1; i++) 
  {  
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  radar_samples_total = 0;
  lidar_samples_total = 0;
  radar_samples_over_NIS_high_percenticle = 0;
  lidar_samples_over_NIS_high_percenticle = 0;
}

UKF::~UKF() {}

void UKF::Initialize(MeasurementPackage measurement_pack)
{

  time_us_ = measurement_pack.timestamp_;
  
  VectorXd p_raw = measurement_pack.raw_measurements_;
  x_.fill(0.0);

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {

    float rho = p_raw(0);
    float phi = p_raw(1);
    float rho_dot = p_raw(2);

    float vx = rho_dot * cos(phi);
    float vy = rho_dot * sin(phi);
    float v_state = sqrt(vx*vx + vy*vy);
    float si = atan2(vx, vy);
    float si_dot = 0.0;

    x_(0) = rho * cos(phi);
    x_(1) = rho * sin(phi);
    x_(2) = v_state;
    x_(3) = si;
    x_(4) = si_dot;


    std::cout << "x,y" << x_(0) << "," << x_(1);

  }
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
  {
    x_(0) = p_raw(0);
    x_(1) = p_raw(1);
  }

  cout << "init was: " << x_(0) << " " << x_(1);
  // done initializing, no need to predict or update
  is_initialized_ = true;
  return;
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
  if(!is_initialized_)
  {
      Initialize(meas_package);
      return;
  }

  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict
  Prediction(dt);

  // Update
  if(meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
  {
    UpdateLidar(meas_package);
  }
  else
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
  // L14 - create augmented sigma points
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  // augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  // augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.fill(0.0);

  // augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  MatrixXd proc_noise = MatrixXd(2,2);
  proc_noise << std_a_ * std_a_, 0,
         0, std_yawdd_ * std_yawdd_;

  // augmented covariance  
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = proc_noise;
  MatrixXd A_aug = P_aug.llt().matrixL();

  // augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_aug_; i++)
  {
      Xsig_aug.col(i+1) = x_aug + (sqrt(lambda_ + n_aug_) * A_aug.col(i));
      Xsig_aug.col(i+1+n_aug_) = x_aug - (sqrt(lambda_ + n_aug_) * A_aug.col(i));
  }

  // Augmented sigma points created... predict sigma points
  for(int i = 0; i < 2*n_aug_+1; i++)
  {
      float px      = Xsig_aug(0, i);
      float py      = Xsig_aug(1, i);
      float v       = Xsig_aug(2, i);
      float u       = Xsig_aug(3, i);
      float u_dot   = Xsig_aug(4, i);
      float va      = Xsig_aug(5, i);
      float v_dot   = Xsig_aug(6, i);
      
      float delta_t_squared = (0.5 * (delta_t * delta_t));
      //std::cout << px << " " << py << " " << v << " " << u << " " << u_dot << " " << va << std::endl; 
      if(fabs(u_dot) < 0.001)
      {
          Xsig_pred_(0,i) = (v * cos(u) * delta_t) + (delta_t_squared * cos(u) * va);
          Xsig_pred_(1,i) = (v * sin(u) * delta_t) + (delta_t_squared * sin(u) * va);

      }
      else
      {
          Xsig_pred_(0,i) = (v/u_dot)*(sin(u + u_dot * delta_t) - sin(u));
          Xsig_pred_(0,i) += (delta_t_squared * cos(u) * va);
          
          Xsig_pred_(1,i) = (v/u_dot)*((-1.0 * cos(u + u_dot * delta_t)) + cos(u));
          Xsig_pred_(1,i) += (delta_t_squared * sin(u) * va);
      }
    
      Xsig_pred_(2,i) = delta_t * va;
      Xsig_pred_(3,i) = (delta_t * u_dot) + delta_t_squared * v_dot;
      Xsig_pred_(4,i) = (delta_t * v_dot);
      
      Xsig_pred_(0,i) += px;
      Xsig_pred_(1,i) += py;
      Xsig_pred_(2,i) += v;
      Xsig_pred_(3,i) += u;
      Xsig_pred_(4,i) += u_dot;
  }

  // Predict Mean and Covariance
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for(int i = 0; i < (2*n_aug_+1); i++)
  {
      //std::cout << "elm: " << std::endl << Xsig_pred(0,i) << std::endl;
      x_pred = x_pred + weights_(i)*Xsig_pred_.col(i);
      
  }


  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  for(int i = 0; i < (2*n_aug_+1); i++)
  {
      VectorXd x_diff_temp = Xsig_pred_.col(i) - x_pred;

      // normalize
      while(x_diff_temp(3) > M_PI)
      {

        x_diff_temp(3) -= 2.0 * M_PI;
      }

      while(x_diff_temp(3) < -M_PI)
      {

        x_diff_temp(3) += 2.0 * M_PI;
      }


      P_pred = P_pred + weights_(i) * (x_diff_temp) * (x_diff_temp.transpose());
    // KDA don't forget to normalize angle difference     
  }

  // Update mean and covariance with predicitions
  x_ = x_pred;
  P_ = P_pred;

  //std::cout << "Prediction" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << P_ << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  const int n_z = 2; // lidar has x, y
  
  //transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;  //x
    Zsig(1,i) = p_y;  //y
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

   //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_lidar_;

  // Measurement Update
  VectorXd z = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);


  // Cross correlation matrix (L30)
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  { 
    VectorXd z_diff = Zsig.col(i) - z_pred;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd K = Tc * S.inverse();
  VectorXd z_delta = z - z_pred;

  // Update state mean and covariance
  x_ = x_ + K * z_delta;
  P_ = P_ - K * S * K.transpose();

  //std::cout << "Update LIDAR" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << P_ << std::endl;

  float NIS_LIDAR = z_delta.transpose() * S.inverse() * z_delta;
  std::cout << "NIS (LIDAR)/expected(95%)" << NIS_LIDAR << "/" << "5.991" << std::endl;
  lidar_samples_total++;

  if(NIS_LIDAR > 5.991)
  {
      lidar_samples_over_NIS_high_percenticle++;
  }
  std::cout << "NIS ratio above 95th percentile " 
            << lidar_samples_over_NIS_high_percenticle 
            << "/" 
            << lidar_samples_total 
            << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  //transform sigma points into measurement space
  const int n_z = 3; // radar has r, phi, r_dot

  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  //2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    
    // watch for divide by zero...
    float denominator = p_x*p_x + p_y*p_y;
    //Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(denominator);   //r_dot
   if(fabs(denominator) < .001)
    {
        Zsig(2,i) = (p_x*v1 + p_y*v2) / .001;
    }
    else
    {
      Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(denominator);   //r_dot
    }
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI){ z_diff(1)-=2.*M_PI;}
    while (z_diff(1)<-M_PI){ z_diff(1)+=2.*M_PI;}

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix

  S = S + R_radar_;

  // Measurement Update
  VectorXd z = meas_package.raw_measurements_;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);


  // Cross correlation matrix (L30)
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  { 
    VectorXd z_diff = Zsig.col(i) - z_pred;

    while (z_diff(1)> M_PI){ z_diff(1)-=2.*M_PI;}
    while (z_diff(1)<-M_PI){ z_diff(1)+=2.*M_PI;}

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    while (x_diff(1)> M_PI){ x_diff(1)-=2.*M_PI;}
    while (x_diff(1)<-M_PI){ x_diff(1)+=2.*M_PI;}

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }



  MatrixXd K = Tc * S.inverse();
  VectorXd z_delta = z - z_pred;

  while (z_delta(1)> M_PI){ z_delta(1)-=2.*M_PI;}
  while (z_delta(1)<-M_PI){ z_delta(1)+=2.*M_PI;}

  // Update state mean and covariance
  x_ = x_ + K * z_delta;

  P_ = P_ - K * S * K.transpose();

  //std::cout << "Update RADAR" << std::endl;
  //std::cout << x_ << std::endl;
  //std::cout << P_ << std::endl;

  float NIS_RADAR = z_delta.transpose() * S.inverse() * z_delta;
  std::cout << "NIS (RADAR)/expected(95%)" << NIS_RADAR << "/" << "7.815" << std::endl;
  radar_samples_total++;

  if(NIS_RADAR > 7.815)
  {
      radar_samples_over_NIS_high_percenticle++;
  }
  std::cout << "NIS (RADAR) ratio above 95th percentile " 
            << radar_samples_over_NIS_high_percenticle 
            << "/" 
            << radar_samples_total 
            << std::endl;
}
