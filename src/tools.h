#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth);

  // helper to generate augmented sigma points
  void GenerateAugmentedSigmaPoints(const VectorXd& state, const MatrixXd& P, const double std_a, const double std_yawdd, const int n_x, const unsigned n_aug, const int lambda, MatrixXd& Xsig_aug);

  // helper to predict sigma points
  void PredictSigmaPoints(const int n_x, const unsigned n_aug, const MatrixXd& Xsig_aug, const double delta_t, MatrixXd& Xsig_pred);

  // helper function to compute mean and covariance of predicted state.
  void PredictNewStateAndStateCovariance(const MatrixXd& Xsig_pred, const VectorXd& weights, VectorXd& state, MatrixXd& P);

  // helper to generate weights to be used for predicting new state and associated covariance matrix
  void GenerateWeights(const unsigned n_aug, const int lambda, VectorXd& weights);

  // helper function to translate predicted sigma points from state space to radar measurement space
  void TranslatePredictedSigmaPointsToRadarMeasurementSpace(const unsigned n_aug, const unsigned nRadarDimension, const MatrixXd& Xsig_pred, MatrixXd& Zsig_pred);


  // helper function to translate predicted sigma points from state space to laser measurement space
  void TranslatePredictedSigmaPointsToLaserMeasurementSpace(const unsigned n_aug, const unsigned nRadarDimension, const MatrixXd& Xsig_pred, MatrixXd& Zsig_pred);

  // helper function to estimate radar measurement based on predicted state
  void EstimateNewRadarMeasurementAndCovarianceMatrix(const double std_radr, const double std_radphi, const double std_radrd, const MatrixXd& predictedSigmaPointsinRadarSpace, const VectorXd& weights, VectorXd& estimatedRadarMeasurement, MatrixXd& S);

  // helper function to estimate laser measurement based on predicted state
  void EstimateNewLaserMeasurementAndCovarianceMatrix(const double std_laspx, const double std_laspy, const MatrixXd& predictedSigmaPointsinLaserSpace, const VectorXd& weights, VectorXd& estimatedLaserMeasurement, MatrixXd& S);

  // helper for compute Kalman gain for radar measurement
  void ComputeKalmanGainForRadar(const MatrixXd& Xsig_pred, const VectorXd& state, const MatrixXd& predictedSigmaPointsinRadarSpace, const VectorXd& estimatedRadarMeasurement, const MatrixXd& S_Radar, const VectorXd& weights, MatrixXd& KalmanGainforRadar);


  // helper for computing Kalman gain for laser case
  void ComputeKalmanGainForLaser(const MatrixXd& Xsig_pred, const VectorXd& state, const MatrixXd& predictedSigmaPointsinLaserSpace, const VectorXd& estimatedLaserMeasurement, const MatrixXd& S_Laser, const VectorXd& weights, MatrixXd& KalmanGainforLaser);

  // helper to compute posteriori state and associated covariance matrix
  void ComputePosterioriStateAndStateCovariance(const MatrixXd& KalmanGainfor, const VectorXd& actualMeasurement, const VectorXd& estimatedMeasurement, const MatrixXd& S, VectorXd& state, MatrixXd& stateCovarianceMatrix);

  //helper to compute NIS for Radar
  void ComputeNISForRadar(const VectorXd& radarMeasurement, const VectorXd& equivalentRadarMeasurement, const MatrixXd& S, double& NISRadar);

  // helper to compute NIS for laser
  void ComputeNISForLaser(const VectorXd& laserMeasurement, const VectorXd& equivalentLaserMeasurement, const MatrixXd& S, double& NISLaser);


private:

  Eigen::VectorXd computeColumnforDeterministicVector(const Eigen::VectorXd& state, const int n_x, const double dt);

  Eigen::VectorXd computeColumnforNonDeterministicVector(const Eigen::VectorXd& state, const int n_x, const double dt);

  double NormalizeYawAngle(double angle);

  VectorXd convertColumnFromStateVectorToRadarMeasurementVector(const unsigned nRadarDimensions, const VectorXd& state);

  VectorXd convertColumnFromStateVectorToLaserMeasurementVector(const unsigned nLaserDimensions, const VectorXd& state);



};

#endif /* TOOLS_H_ */
