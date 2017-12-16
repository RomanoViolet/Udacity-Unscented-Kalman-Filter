#include <iostream>
#include "tools.h"
#include <cmath>
#include <functional>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
}

Tools::~Tools() {
}

/* -----------------------------------------------------
 * ------------------  PUBLIC METHODS ------------------
 * -----------------------------------------------------
 */


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {
	/**
	 TODO:
	 * Calculate the RMSE here.
	 */

	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	eigen_assert(
			(estimations.size() != 0)
					&& (estimations.size() == ground_truth.size())
					&& "size of vectors do not match.");

	//accumulate squared residuals
	VectorXd totalError(4);
	totalError << 0, 0, 0, 0;

	for (unsigned i = 0; i < estimations.size(); ++i) {
		VectorXd currentEstimate = estimations[i];
		VectorXd currentGroundTruth = ground_truth[i];
		VectorXd currentError = currentEstimate - currentGroundTruth;
		totalError = totalError
				+ static_cast<VectorXd>(currentError.array().pow(2.0));
	}

	//calculate the mean
	VectorXd meanError = totalError.array() / estimations.size();

	//calculate the squared root
	rmse = meanError.array().pow(0.5);

	//return the result
	return rmse;

}/*Tools::CalculateRMSE*/


void Tools::GenerateAugmentedSigmaPoints(const VectorXd& state, const MatrixXd& P, const double std_a, const double std_yawdd, const int n_x, const unsigned n_aug, const int lambda, MatrixXd& Xsig_aug)
{

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug);

	// Q vector
	MatrixXd Q = MatrixXd(n_aug - n_x, n_aug - n_x);
	Q << std::pow(std_a, 2.0), 0, 0, std::pow(std_yawdd, 2.0);

	//create augmented state covariance -- zero initialized.
	MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);

	// fill in the augmented matrix
	// upper left half
	P_aug.block(0, 0, n_x, n_x) = P;

	// lower right half
	P_aug.block(n_x, n_x, n_aug - n_x, n_aug - n_x) = Q;



	//calculate square root of P
	MatrixXd A = P_aug.llt().matrixL();


	//create augmented mean state
	VectorXd means = state.rowwise().mean();


	//the mean states of the process noise
	Eigen::VectorXd meansProcessNoise = Eigen::VectorXd(2);

	// noise parameters have zero mean
	meansProcessNoise << 0, 0;

	// fill the first (i.e., means) column of the augmented sigma point matrix
	Xsig_aug.col(0) << means, meansProcessNoise;


	VectorXd adder = VectorXd::Zero(n_x);
	adder = Xsig_aug.col(0);


	// left half of augmented matrix (right of means column)
	MatrixXd leftHalf = (std::pow((lambda + n_aug), 0.5) * A).array().colwise()
			+ Xsig_aug.col(0).array();


	// right half of augmented matrix where Xsig_aug.col(0) is means vector
	MatrixXd rightHalf =
			((-1 * std::pow((lambda + n_aug), 0.5)) * A).array().colwise()
					+ Xsig_aug.col(0).array();


	// populate the rest of augmented sigma point matrix
	Xsig_aug.block(0, 1, n_aug, n_aug) = leftHalf;
	Xsig_aug.block(0, n_aug + 1, n_aug, n_aug) = rightHalf;


}/*Tools::GenerateAugmentedSigmaPoints*/




void Tools::PredictSigmaPoints(const int n_x, const unsigned n_aug, const MatrixXd& Xsig_aug, const double delta_t, MatrixXd& Xsig_pred)
{

	  Xsig_pred = MatrixXd::Zero(n_x, 2 * n_aug + 1);
	  VectorXd columnDeterministicPart = VectorXd::Zero(n_x);
	  VectorXd columnNonDeterministicPart = VectorXd::Zero(n_x);
	  for (size_t columnCounter = 0; columnCounter < Xsig_aug.cols(); columnCounter++)
	  {
	    columnDeterministicPart = this->computeColumnforDeterministicVector(Xsig_aug.col(columnCounter), n_x, delta_t);
	    columnNonDeterministicPart = this->computeColumnforNonDeterministicVector(Xsig_aug.col(columnCounter), n_x, delta_t);

	    Xsig_pred.col(columnCounter) = columnDeterministicPart + columnNonDeterministicPart + Xsig_aug.col(columnCounter).segment(0,n_x);

	  }

}/*Tools::PredictSigmaPoints*/



void Tools::GenerateWeights(const unsigned n_aug, const int lambda, VectorXd& weights)
{
	//special case for index = 0
	unsigned n_sig = 2 * n_aug + 1;

	//zero the weights
	weights.Zero(weights.cols());

	//weights(0) = (lambda*1.0)/(lambda+n_sig);
	weights(0) = (lambda * 1.0) / (lambda + n_aug);

	// for elements 1:(2*n_aug+1)
	weights.tail(n_sig - 1).setConstant((0.5) / (lambda + n_aug));

}/*Tools::GenerateWeights*/


void Tools::PredictNewStateAndStateCovariance(const MatrixXd& Xsig_pred, const VectorXd& weights, VectorXd& state, MatrixXd& P)
{

	// predict the new state
	state = Xsig_pred * weights;


	// compute the predicted state covariance matrix, P
	MatrixXd P_1 = Xsig_pred.array().colwise() - state.array();


	// normalize the yaw angle by providing a functor to non-static NormalizeYawAngle function
	P_1.row(3) = P_1.row(3).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	// multiply each column of P_1 with the corresponding column of weights.transpose().
	// That is, weights.transpose() is a 1x15 vector
	MatrixXd P_2 = P_1.array().rowwise() * weights.transpose().array();


	// final state covariance matrix
	P = P_2 * P_1.transpose();

}/* Tools::PredictNewStateAndStateCovariance */


void Tools::TranslatePredictedSigmaPointsToRadarMeasurementSpace(const unsigned n_aug, const unsigned nRadarDimension, const MatrixXd& Xsig_pred, MatrixXd& Zsig_pred)
{

	//Zsig_pred = Xsig_pred;
	// Convert state-state to measurement space, column-by-column (i.e., create Z (upper case "z")
	for (size_t column = 0; column < 2 * n_aug + 1; column++) {
		// convert each column from predicted state-space (Xsig_pred) to measurement space ("Z", or here, Zsig) individually.
		Zsig_pred.col(column) =
				this->convertColumnFromStateVectorToRadarMeasurementVector(
						nRadarDimension, Xsig_pred.col(column));
	}

}/*Tools::TranslatePredictedSigmaPointsToRadarMeasurementSpace*/



void Tools::TranslatePredictedSigmaPointsToLaserMeasurementSpace(const unsigned n_aug, const unsigned nLaserDimension, const MatrixXd& Xsig_pred, MatrixXd& Zsig_pred)
{
	// Convert state-state to measurement space, column-by-column (i.e., create Z (upper case "z")
	for (size_t column = 0; column < 2 * n_aug + 1; column++) {
		// convert each column from predicted state-space (Xsig_pred) to measurement space ("Z", or here, Zsig) individually.
		Zsig_pred.col(column) =
				this->convertColumnFromStateVectorToLaserMeasurementVector(
						nLaserDimension, Xsig_pred.col(column));
	}
}/*Tools::TranslatePredictedSigmaPointsToLaserMeasurementSpace*/


void Tools::EstimateNewRadarMeasurementAndCovarianceMatrix(const double std_radr, const double std_radphi, const double std_radrd, const MatrixXd& predictedSigmaPointsinRadarSpace, const VectorXd& weights, VectorXd& estimatedRadarMeasurement, MatrixXd& S)
{
	// initialization
	MatrixXd Zsig_Weighted = MatrixXd::Zero(
			predictedSigmaPointsinRadarSpace.rows(),
			predictedSigmaPointsinRadarSpace.cols());

	Zsig_Weighted = predictedSigmaPointsinRadarSpace.array().rowwise()
			* weights.transpose().array();

	// normalize
	Zsig_Weighted.row(1) = Zsig_Weighted.row(1).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));


	// declare a helper vector which will be used to sum of all elements of each row in Zsig_weighted to yield one element
	VectorXd ones = VectorXd::Ones(predictedSigmaPointsinRadarSpace.cols());

	// compute the final predicted measurement based on current state: 3x1
	estimatedRadarMeasurement = Zsig_Weighted * ones;

	// normalize again
	estimatedRadarMeasurement.row(1) = estimatedRadarMeasurement.row(1).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));


	MatrixXd S_1 = predictedSigmaPointsinRadarSpace.array().colwise()
			- estimatedRadarMeasurement.array();

	/*Normalize angles*/
	S_1.row(1) = S_1.row(1).unaryExpr(
			std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	//std::cout << "S_1" << std::endl << S_1 << std::endl;

	//Multiply weights columnwise (i.e., each column gets the same weight) from the corresponding weight matrix
	MatrixXd S_2 = S_1.array().rowwise() * weights.transpose().array();

	// normalize angles
	S_2.row(1) = S_2.row(1).unaryExpr(
				std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	// And then multiply it with the transpose. S_3: 3x3
	MatrixXd S_3 = S_2 * S_1.transpose();
	//std::cout << "S_3" << std::endl << S_3 << std::endl << std::endl;

	// Add noise (i.e., non-deterministic) component
	MatrixXd R = MatrixXd::Zero(predictedSigmaPointsinRadarSpace.rows(),
			predictedSigmaPointsinRadarSpace.rows());
	R(0, 0) = std::pow(std_radr, 2.0);
	R(1, 1) = std::pow(std_radphi, 2.0);
	R(2, 2) = std::pow(std_radrd, 2.0);

	// Final measurement covariance matrix S
	S = S_3 + R;

}/*Tools::EstimateNewRadarMeasurementAndCovarianceMatrix*/


void Tools::EstimateNewLaserMeasurementAndCovarianceMatrix(const double std_laspx, const double std_laspy, const MatrixXd& predictedSigmaPointsinLaserSpace, const VectorXd& weights, VectorXd& estimatedLaserMeasurement, MatrixXd& S)
{

	// initialization
	MatrixXd Zsig_Weighted = MatrixXd::Zero(
			predictedSigmaPointsinLaserSpace.rows(),
			predictedSigmaPointsinLaserSpace.cols());


	Zsig_Weighted = predictedSigmaPointsinLaserSpace.array().rowwise()
			* weights.transpose().array();

	// declare a helper vector which will be used to sum of all elements of each row in Zsig_weighted to yield one element
	VectorXd ones = VectorXd::Ones(predictedSigmaPointsinLaserSpace.cols());

	// compute the final predicted measurement based on current state: 2x1
	estimatedLaserMeasurement = Zsig_Weighted * ones;


	MatrixXd S_1 = predictedSigmaPointsinLaserSpace.array().colwise()
			- estimatedLaserMeasurement.array();



	//Multiply weights columnwise (i.e., each column gets the same weight) from the corresponding weight matrix
	MatrixXd S_2 = S_1.array().rowwise() * weights.transpose().array();


	// And then multiply it with the transpose. S_3: 3x3
	MatrixXd S_3 = S_2 * S_1.transpose();


	// Add noise (i.e., non-deterministic) component
	MatrixXd R = MatrixXd::Zero(predictedSigmaPointsinLaserSpace.rows(),
			predictedSigmaPointsinLaserSpace.rows());



	R(0, 0) = std::pow(std_laspx, 2.0);
	R(1, 1) = std::pow(std_laspy, 2.0);

	// Final measurement covariance matrix S
	S = S_3 + R;


}/*Tools::EstimateNewLaserMeasurementAndCovarianceMatrix*/



void Tools::ComputeKalmanGainForRadar(const MatrixXd& Xsig_pred, const VectorXd& state, const MatrixXd& predictedSigmaPointsinRadarSpace, const VectorXd& estimatedRadarMeasurement, const MatrixXd& S_Radar, const VectorXd& weights, MatrixXd& KalmanGainforRadar)
{

	MatrixXd Tc = MatrixXd(state.cols(),
			predictedSigmaPointsinRadarSpace.rows());

	// First compute w(X-x)s
	MatrixXd T_1 = Xsig_pred.array().colwise() - state.array();
	//std::cout << "P_1: " << std::endl << P_1 << std::endl;

	T_1.row(3) = T_1.row(3).unaryExpr(
			std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	// multiply each column of T_1 with the corresponding column of weights.transpose().
	// That is, weights.transpose() is a 1x15 vector
	MatrixXd T_2 = T_1.array().rowwise() * weights.transpose().array();

	// Compute Z-z
	MatrixXd T_3 = predictedSigmaPointsinRadarSpace.array().colwise()
			- estimatedRadarMeasurement.array();

	// normalize
	T_3.row(1) = T_3.row(1).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	// Compute the product of T_2 and T_3
	Tc = T_2 * T_3.transpose();

	// The Kalman gain
	KalmanGainforRadar = Tc * S_Radar.inverse();


}/*Tools::ComputeKalmanGainForRadar*/




void Tools::ComputeKalmanGainForLaser(const MatrixXd& Xsig_pred, const VectorXd& state, const MatrixXd& predictedSigmaPointsinLaserSpace, const VectorXd& estimatedLaserMeasurement, const MatrixXd& S_Laser, const VectorXd& weights, MatrixXd& KalmanGainforLaser)
{

	MatrixXd Tc = MatrixXd::Zero(state.rows(),
			predictedSigmaPointsinLaserSpace.rows());

	// First compute w(X-x)s
	MatrixXd T_1 = Xsig_pred.array().colwise() - state.array();
	T_1.row(3) = T_1.row(3).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));
	// normalize


	// multiply each column of T_1 with the corresponding column of weights.transpose().
	// That is, weights.transpose() is a 1x15 vector
	MatrixXd T_2 = T_1.array().rowwise() * weights.transpose().array();

	// normalize again
	T_2.row(3) = T_2.row(3).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));


	// Compute Z-z
	MatrixXd T_3 = predictedSigmaPointsinLaserSpace.array().colwise()
			- estimatedLaserMeasurement.array();


	// Compute the product of T_2 and T_3
	Tc = T_2 * T_3.transpose();


	// The Kalman gain
	KalmanGainforLaser = Tc * S_Laser.inverse();

}/*Tools::ComputeKalmanGainForLaser*/


void Tools::ComputePosterioriStateAndStateCovariance(const MatrixXd& KalmanGain, const VectorXd& actualMeasurement, const VectorXd& estimatedMeasurement, const MatrixXd& S, VectorXd& state, MatrixXd& stateCovarianceMatrix)
{
	// Updated state estimate taking current measurement into account (posterior)
	// x: 5x1, K: 5x3, z(current measurement from the sensor): 3x1, z_pred(predicted measurement): 3x1
	state = state
			+ KalmanGain
					* (actualMeasurement - estimatedMeasurement);

	// normalize
	//state.row(3) = state.row(3).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	// Update state covariance matrix. P: 5x5, K: 5x3, S: 3x3, K': 3x5
	stateCovarianceMatrix = stateCovarianceMatrix
			- KalmanGain * S * KalmanGain.transpose();

}/*Tools::ComputePosterioriStateAndStateCovariance*/



void Tools::ComputeNISForRadar(const VectorXd& radarMeasurement, const VectorXd& equivalentRadarMeasurement, const MatrixXd& S, double& NISRadar)
{
	//std::cout << "radarMeasurement" << std::endl << radarMeasurement << std::endl;
	//std::cout << "equivalentRadarMeasurement" << std::endl << equivalentRadarMeasurement << std::endl;
	//std::cout << "S" << std::endl << S << std::endl;
	VectorXd error = radarMeasurement - equivalentRadarMeasurement;
	/*Normalize yaw angle*/
	error.row(1) = error.row(1).unaryExpr(std::bind(&Tools::NormalizeYawAngle, this, std::placeholders::_1));

	NISRadar = error.transpose() * S.inverse() * error;

}/*Tools::ComputeNISForRadar*/



void Tools::ComputeNISForLaser(const VectorXd& laserMeasurement, const VectorXd& equivalentLaserMeasurement, const MatrixXd& S, double& NISLaser)
{
	//std::cout << "laserMeasurement" << std::endl << laserMeasurement << std::endl;
	//std::cout << "equivalentLaserMeasurement" << std::endl << equivalentLaserMeasurement << std::endl;
	//std::cout << "S" << std::endl << S << std::endl;

	NISLaser = (laserMeasurement - equivalentLaserMeasurement).transpose() * S.inverse() * (laserMeasurement - equivalentLaserMeasurement);

}/*Tools::ComputeNISForLaser*/




/* -----------------------------------------------------
 * --------------- END OF PUBLIC METHODS ---------------
 * -----------------------------------------------------
 */







/* -----------------------------------------------------
 * ------------------  PRIVATE METHODS ------------------
 * -----------------------------------------------------
 */


double Tools::NormalizeYawAngle(double angle) {


	while (angle>  M_PI) angle-=(2.0 * M_PI);
	while (angle< -M_PI) angle+= (2.0 * M_PI);
	return (angle);

	if (std::abs(angle) / (M_PI) > 1) {
		//https://stackoverflow.com/questions/24234609/standard-way-to-normalize-an-angle-to-%CF%80-radians-in-java
		//theta - TWO_PI * Math.floor((theta + Math.PI) / TWO_PI)

		//version 1
		//return (angle - (std::floor((angle + M_PI) / (2 * M_PI)) * 2 * M_PI));

		//version 2
		return (std::atan2(std::sin(angle), std::cos(angle)));
	} else {
		return (angle);
	}

}/* Tools::NormalizeYawAngle */


Eigen::VectorXd Tools::computeColumnforNonDeterministicVector(const Eigen::VectorXd& state, const int n_x, const double dt)
{
  VectorXd columnWiseResult = VectorXd::Zero(n_x);
  double yawAngle = state(3);
  double accelerationNoise = state(5);
  double yawRateNoise = state(6);
  columnWiseResult(0) = 0.5*(std::pow(dt, 2.0))*std::cos(yawAngle)*accelerationNoise;
  columnWiseResult(1) = 0.5*(std::pow(dt, 2.0))*std::sin(yawAngle)*accelerationNoise;
  columnWiseResult(2) = dt*accelerationNoise;
  columnWiseResult(3) = this->NormalizeYawAngle(0.5*(std::pow(dt, 2.0))*yawRateNoise);
  columnWiseResult(4) = dt*yawRateNoise;

  return columnWiseResult;

}/* Tools::computeColumnforNonDeterministicVector */




Eigen::VectorXd Tools::computeColumnforDeterministicVector(const Eigen::VectorXd& state, const int n_x, const double dt)
{
  VectorXd columnWiseResult = VectorXd::Zero(n_x);

  double velocity = state(2);
  double yawAngle = state(3);
  double yawRate = state(4);
  if(std::fabs(yawRate)>0.001)
  {
    // yaw rate is non-zero
      columnWiseResult(0) = (velocity/yawRate)*(std::sin(yawAngle + yawRate*dt) - std::sin(yawAngle));
      columnWiseResult(1) = (velocity/yawRate)*(-std::cos(yawAngle + yawRate*dt) + std::cos(yawAngle));
      columnWiseResult(2) = 0;
      columnWiseResult(3) = this->NormalizeYawAngle(yawRate*dt);
      columnWiseResult(4) = 0;

  }
  else
  {
      columnWiseResult(0) = (velocity)*std::cos(yawAngle)*dt;
      columnWiseResult(1) = (velocity)*std::sin(yawAngle)*dt;
      columnWiseResult(2) = 0;
      columnWiseResult(3) = this->NormalizeYawAngle(yawRate*dt);
      columnWiseResult(4) = 0;
  }


  /* Normalize Yaw Angle*/

  return columnWiseResult;

}/* Tools::computeColumnforDeterministicVector */



VectorXd Tools::convertColumnFromStateVectorToLaserMeasurementVector(const unsigned nLaserDimensions, const VectorXd& state)
{
	double px = state(0);
	double py = state(1);

	VectorXd result = VectorXd::Zero(nLaserDimensions);
	result << px, py;

	return (result);
}/* Tools::convertColumnFromStateVectorToLaserMeasurementVector */





VectorXd Tools::convertColumnFromStateVectorToRadarMeasurementVector(const unsigned nRadarDimensions, const VectorXd& state)
{
	double px = state(0);
	double py = state(1);
	double velocity = state(2);
	double yawAngle = state(3);

	// radial distance from radar
	double rho = std::pow((std::pow(px, 2.0) + std::pow(py, 2.0)), 0.5);

	// orientation from radar
	double phi = this->NormalizeYawAngle(std::atan2(py, px));

	// radial velocity: rhodot
	double rhodot = ((px * std::cos(yawAngle) * velocity)
			+ (py * std::sin(yawAngle) * velocity))
			/ std::pow((std::pow(px, 2.0) + std::pow(py, 2.0)), 0.5);

	VectorXd result = VectorXd::Zero(nRadarDimensions);
	result << rho, phi, rhodot;

	return (result);
}/* Tools::convertColumnFromStateVectorToRadarMeasurementVector */




/* -----------------------------------------------------
 * --------------- END OF PRIVATE METHODS ---------------
 * -----------------------------------------------------
 */

