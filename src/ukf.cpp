#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <assert.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// if this is false, laser measurements will be ignored (except during init)
	use_laser_ = true;

	// if this is false, radar measurements will be ignored (except during init)
	use_radar_ = true;

	// Process noise standard deviation longitudinal acceleration in m/s^2
	/*
	 */
	std_a_ = 0.6;


	// Process noise standard deviation yaw acceleration in rad/s^2
	/*
	 * The picture from lesson indicates the value to be 0.5
	 */
	std_yawdd_ = 0.5;
	//std_yawdd_ = 0.1;

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

	// the dimension of state vector
	n_x_ = 5;

	// dimension of the augmented state vector
	n_aug_ = 7;

	// initial state vector
	x_ = VectorXd(n_x_);

	// initial covariance matrix
	P_ = MatrixXd::Zero(n_x_, n_x_);
	
	// initialize P_
	P_(0, 0) = 1.0;
	P_(1, 1) = 1.0;
	P_(2, 2) = 1.0;
	P_(3, 3) = 1.0;
	P_(4, 4) = 1.0;

	// initialize augmented sigma point matrix
	Xsig_aug_ = MatrixXd::Zero(n_aug_, 2 * n_aug_ + 1);


	//define spreading parameter
	lambda_ = static_cast<double>(static_cast<signed>(3)-static_cast<signed>(n_aug_));

	// initialize weights vector
	weights_ = VectorXd(2 * n_aug_ + 1);

	// initialize time_us_ so that compiler stops complaining
	time_us_ = 0.0F;

	// initialize dt so that compiler stops complaining
	dt = 0.0F;

	// initialize the vector to hold radar measurements
	actualRadarMeasurement = VectorXd::Zero(nRadarDimension);

	// initialize the vector to hold laser measurements
	actualLaserMeasurement = VectorXd::Zero(nLaserDimension);

	// UKF is set to initially not initialized
	is_initialized_ = false;

	// initialize the matrix to hold predicted sigma points translated from state space to radar space
	predictedSigmaPointsinRadarSpace = MatrixXd::Zero(nRadarDimension, 2*n_aug_ + 1);

	// initialize the matrix to hold predicted sigma points translated from state space to laser space
	predictedSigmaPointsinLaserSpace = MatrixXd::Zero(nLaserDimension, 2*n_aug_ + 1);

	// initialize a matrix to hold predicted sigma points
	Xsig_pred_ = MatrixXd::Zero(n_x_, 2 * n_aug_ + 1);

	// initialize a vector to hold estimated radar measurement based on predicted state
	estimatedRadarMeasurement = VectorXd::Zero(nRadarDimension);

	// initialize a vector to hold estimated radar measurement based on predicted state
	estimatedLaserMeasurement = VectorXd::Zero(nLaserDimension);

	// initialize measurement covariance matrix for radar
	S_Radar = MatrixXd::Zero(nRadarDimension, nRadarDimension);

	// initialize measurement covariance matrix for laser
	S_Laser = MatrixXd::Zero(nLaserDimension, nLaserDimension);

	// intialize the matrix to hold kalman gain for radar
	KalmanGainforRadar = MatrixXd::Zero(n_x_, nRadarDimension);

	// intialize the matrix to hold kalman gain for laser
	KalmanGainforLaser = MatrixXd::Zero(n_x_, nLaserDimension);

	// initialize NIS_radar_ so that compiler stops complaining
	NIS_radar_ = 0.0;

	// initialize NIS_laser_ so that compiler stops complaining
	NIS_laser_ = 0.0;

	/**
	 TODO:

	 Complete the initialization. See ukf.h for other member properties.

	 Hint: one or more values initialized above might be wildly off...
	 */
}

UKF::~UKF() {
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_package) {
	/**
	 TODO:

	 Complete this function! Make sure you switch between lidar and radar
	 measurements.
	 */

	/*
	 * Step 1: Initialize state and precision matrices
	 */
	if (!is_initialized_) {

		// the first measurement will be used to update the vectors

		// variable declarations
		double rho;
		double theta;
		//double dtrho;
		double px;
		double py;
		double velocity;
		double yawAngle;
		double yawRate;

		if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {

			/**
			 Convert radar from polar to cartesian coordinates and initialize state.
			 */

			rho = measurement_package.raw_measurements_[0];
			theta = measurement_package.raw_measurements_[1];
			//dtrho = measurement_package.raw_measurements_[2];

			// The inferred px. The supplied angle is in radians.
			px = rho * cos(theta);

			// The inferred py. The supplied angle is in radians.
			py = rho * sin(theta);

			/* The inferred velocity from the object's frame is not directly measured by the radar.
			 *  The radar measures the velocity along the radial distance to the object
			 *  The velocity in the frame of the object is initialized to be double (random number > 1.0) to what radar measures
			 *  to account for the measurement setup, and assuming that the object is indeed moving.
			 */
			velocity = 0.0;

			// The yaw angle: knowledge requires integration over time, so we initially set to zero.
			yawAngle = 0.0;

			// The yaw rate: tje vehicle is not moving at the beginning, so we set to zero.
			yawRate = 0.0;

		} else if (measurement_package.sensor_type_
				== MeasurementPackage::LASER) {

			// px is directly measured by the laser
			px = measurement_package.raw_measurements_[0];

			// py is directly measured by the laser
			py = measurement_package.raw_measurements_[1];

			// laser cannot measure velocity directly: initialize to zero assuming stationary car at the beginning.
			velocity = 0.0;

			// laser cannot measure yaw angle directly: initialize to zero assuming stationary car at the beginning.
			yawAngle = 0.0;

			// laser cannot measure yaw rate directly: initialize to zero assuming stationary car at the beginning.
			yawRate = 0.0;

			// initialize the state vector
			x_ << px, py, velocity, yawAngle, yawRate;



		}/*else if (measurement_package.sensor_type_ == MeasurementPackage::LASER)*/


		// update the time stamp
		time_us_ = measurement_package.timestamp_;

		// set the initialize flag to true
		is_initialized_ = true;

		// do not process further in the initialization round
		return;

	}/*if (!is_initialized_)*/


	// predict sigma points. The result in held in Xsig_pred_
	dt = static_cast<double>((measurement_package.timestamp_ - time_us_) / 1E6);

	// update time stamp
	time_us_ = measurement_package.timestamp_;

	// generate augmented sigma points
	/*
	 * inputs:
	 * 			state vector: x_
	 * 			state covariance matrix: P_
	 * 			acceleration noise (standard deviation): std_a
	 * 			yaw rate noise (standard deviation): std_yawdd_
	 * 			number of dimensions in the state vector (aka rows): n_x_
	 * 			number of dimensions of the augmented state vector (aka rows): n_aug_
	 * 			spreading paramters: lambda_
	 * outputs:
	 *			Augmented sigma points: Xsig_aug_
	 */
	tools.GenerateAugmentedSigmaPoints(x_, P_, std_a_, std_yawdd_,
			n_x_, n_aug_, lambda_, Xsig_aug_);



  // *See detail discussion* @ https://discussions.udacity.com/t/numerical-instability-of-the-implementation/230449/55
	while (dt > 0.1)
	{
		const double constantDt = 0.05;
		this->Prediction(constantDt);
		dt = dt - constantDt;
	}


	this->Prediction(dt);

	// update are sensitive to the sensor from which measurements have been obtained.
	if (measurement_package.sensor_type_ == MeasurementPackage::RADAR) {
		this->UpdateRadar(measurement_package);
	} else {
		this->UpdateLidar(measurement_package);
	}

}/*ProcessMeasurement(MeasurementPackage measurement_package)*/

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	// prediction of sigma points
	/*
	 * inputs:
	 * 			size of state vector: n_x_
	 * 			size of augmented vector: n_aug_
	 * 			augmented matrix: Xsig_aug_
	 * 			time in seconds since last call of UKF: delta_t
	 * outputs:
	 * 			predicted sigma point matrix: Xsig_pred_
	 */
	tools.PredictSigmaPoints(n_x, n_aug_, Xsig_aug_, delta_t, Xsig_pred_);



	// generate weights to be used for computing predicted state and associated
	/*
	 * inputs:
	 * 			size of augmented vector: n_aug_
	 * 			spreading parameter: lambda_
	 * outputs:
	 * 			vector of weights to go along with predicted sigma points: weights_
	 */
	tools.GenerateWeights(n_aug_, lambda_, weights_);



	// predict new state and associated covariance matrix
	/*
	 * inputs:
	 *			predicted sigma points: Xsig_pred_
	 *			vector of weights to go along with predicted sigma points: weights_
	 * outputs:
	 * 			vector of predicted state: x_
	 * 			associated state covariance matrix: P_
	 */
	tools.PredictNewStateAndStateCovariance(Xsig_pred_, weights_, x_, P_);

}




/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

	// clean the matrix
	predictedSigmaPointsinLaserSpace.setZero();

	// translate predicted sigma points into measurement space
	/*
	 * inputs
	 * 			size of augmented vector: n_aug_
	 * 			size of augmented vector: n_aug_
	 * 			predicted sigma point matrix: Xsig_pred_
	 * outputs
	 * 			matrix to hold sigma points translated from state space to radar measurement space: predictedSigmaPointsinRadarSpace
	 */
	tools.TranslatePredictedSigmaPointsToLaserMeasurementSpace(n_aug_, nLaserDimension, Xsig_pred_, predictedSigmaPointsinLaserSpace);



	// clean up the vector
	estimatedLaserMeasurement.setZero();
	// cleanup the measurement covariance matrix
	S_Laser.setZero();

	// estimate radar measurement based on predicted state space, and the corresponding measurement covariance matrix
	/*
	 * inputs
	 * 			standard deviation in the px component of measurement: std_laspx_
	 *  		standard deviation in the py component of measurement: std_laspy_
	 * 			matrix with translated from state space to laser measurement space: predictedSigmaPointsinLaserSpace
	 * 			vector of weights to go along with predicted sigma points: weights_
	 * outputs
	 * 			vector with estimated laser measurement based on predicted state: estimatedLaserMeasurement
	 * 			matrix with associated measurement covariances: S_Laser
	 */
	tools.EstimateNewLaserMeasurementAndCovarianceMatrix(std_laspx_, std_laspy_, predictedSigmaPointsinLaserSpace, weights_, estimatedLaserMeasurement, S_Laser);



	// cleanup previous kalman gain
	KalmanGainforLaser.setZero();
	// compute kalman gain
	/*
	 * inputs:
	 * 			Matrix consisting of predicted sigma points in state space: Xsig_pred_
	 * 			Predicted state: x_
	 * 			Matrix consisting of predicted sigma points translated in laser space: predictedSigmaPointsinLaserSpace
	 * 			Vector with estimated measurement based on predicted state: estimatedLaserMeasurement
	 * 			Estimated Measurement covariance matrix: S_Laser
	 * 			vector of weights to go along with predicted sigma points: weights_
	 * outputs:
	 * 	 		Kalman Gain in the case of laser measurement: KalmanGainforLaser
	 */
	tools.ComputeKalmanGainForLaser(Xsig_pred_, x_, predictedSigmaPointsinLaserSpace, estimatedLaserMeasurement, S_Laser, weights_, KalmanGainforLaser);




	// fill in the laser measurements into the vector
	actualLaserMeasurement << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

	// compute posteriori state and convariance estimate
	/*
	 * inputs:
	 * 			Kalman Gain in the case of laser measurement: KalmanGainforLaser
	 * 			Actual incoming measurement for the laser: actualLaserMeasurement
	 * 			Vector with estimated measurement based on predicted state: estimatedLaserMeasurement
	 * 			Estimated Measurement covariance matrix: S_Laser
	 * 	inputs and outputs (in-situ modified variables)
	 * 			updated state: x_
	 * 			updated state covariance matrix: P_
	 */
	tools.ComputePosterioriStateAndStateCovariance(KalmanGainforLaser, actualLaserMeasurement, estimatedLaserMeasurement, S_Laser, x_, P_);



	// compute NIS values for laser
		/*
		 * inputs:
		 * 			Actual incoming measurement for the radar: actualLaserMeasurement
		 * 			Vector with estimated measurement based on predicted state: estimatedLaserMeasurement
		 * 			Estimated Measurement covariance matrix: S_Laser
		 * outputs
		 * 			NIS for Radar: NIS_laser_
		 */
	tools.ComputeNISForLaser(actualLaserMeasurement, estimatedLaserMeasurement, S_Laser, NIS_laser_);

}/*UKF::UpdateLidar*/





/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

	// clean the matrix
	predictedSigmaPointsinRadarSpace.setZero();

	// translate predicted sigma points into measurement space
	/*
	 * inputs
	 * 			size of augmented vector: n_aug_
	 * 			size of augmented vector: n_aug_
	 * 			predicted sigma point matrix: Xsig_pred_
	 * outputs
	 * 			matrix to hold sigma points translated from state space to radar measurement space: predictedSigmaPointsinRadarSpace
	 */
	tools.TranslatePredictedSigmaPointsToRadarMeasurementSpace(n_aug_, nRadarDimension, Xsig_pred_, predictedSigmaPointsinRadarSpace);



	// clean up the vector
	estimatedRadarMeasurement.setZero();
	// cleanup the measurement covariance matrix
	S_Radar.setZero();

	// estimate radar measurement based on predicted state space, and the corresponding measurement covariance matrix
	/*
	 * inputs
	 * 			standard deviation in the rho component of measurement: std_radr_
	 *  		standard deviation in the phi component of measurement: std_radphi_
	 *  		standard deviation in the rhodot component of measurement: std_radrd_
	 * 			matrix with translated from state space to radar measurement space: predictedSigmaPointsinRadarSpace
	 * 			vector of weights to go along with predicted sigma points: weights_
	 * outputs
	 * 			vector with estimated radar measurement based on predicted state: estimatedRadarMeasurement
	 * 			matrix with associated measurement covariances: S_Radar
	 */
	tools.EstimateNewRadarMeasurementAndCovarianceMatrix(std_radr_, std_radphi_, std_radrd_, predictedSigmaPointsinRadarSpace, weights_, estimatedRadarMeasurement, S_Radar);



	// cleanup previous kalman gain
	KalmanGainforRadar.setZero();
	// compute kalman gain
	/*
	 * inputs:
	 * 			Matrix consisting of predicted sigma points in state space: Xsig_pred_
	 * 			Predicted state: x_
	 * 			Matrix consisting of predicted sigma points translated in radar space: predictedSigmaPointsinRadarSpace
	 * 			Vector with estimated measurement based on predicted state: estimatedRadarMeasurement
	 * 			Estimated Measurement covariance matrix: S_Radar
	 * 			vector of weights to go along with predicted sigma points: weights_
	 * outputs:
	 * 	 		Kalman Gain in the case of radar measurement: KalmanGainforRadar
	 */
	tools.ComputeKalmanGainForRadar(Xsig_pred_, x_, predictedSigmaPointsinRadarSpace, estimatedRadarMeasurement, S_Radar, weights_, KalmanGainforRadar);



	// fill in the radar measurements into the vector
	actualRadarMeasurement << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];



	// compute posteriori state and convariance estimate
	/*
	 * inputs:
	 * 			Kalman Gain in the case of radar measurement: KalmanGainforRadar
	 * 			Actual incoming measurement for the radar: actualRadarMeasurement
	 * 			Vector with estimated measurement based on predicted state: estimatedRadarMeasurement
	 * 			Estimated Measurement covariance matrix: S_Radar
	 * 	inputs and outputs (in-situ modified variables)
	 * 			updated state: x_
	 * 			updated state covariance matrix: P_
	 */
	tools.ComputePosterioriStateAndStateCovariance(KalmanGainforRadar, actualRadarMeasurement, estimatedRadarMeasurement, S_Radar, x_, P_);



	// compute NIS values for radar
	/*
	 * inputs:
	 * 			Actual incoming measurement for the radar: actualRadarMeasurement
	 * 			Vector with estimated measurement based on predicted state: estimatedRadarMeasurement
	 * 			Estimated Measurement covariance matrix: S_Radar
	 * outputs
	 * 			NIS for Radar: NIS_radar_
	 */
	tools.ComputeNISForRadar(actualRadarMeasurement, estimatedRadarMeasurement, S_Radar, NIS_radar_);


}
