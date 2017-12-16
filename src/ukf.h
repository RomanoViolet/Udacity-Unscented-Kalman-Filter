#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

	///* initially set to false, set to true in first call of ProcessMeasurement
	bool is_initialized_;

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

	///* augmented sigma points matrix
	MatrixXd Xsig_aug_;

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
	double std_radrd_;

	///* Weights of sigma points
	VectorXd weights_;

	///* State dimension
	int n_x_;

	///* Augmented state dimension
	unsigned n_aug_;

	///* Sigma point spreading parameter
	double lambda_;

	// NIS for laser
	double NIS_laser_;

	// NIS for radar
	double NIS_radar_;

	//set state dimension
	const unsigned n_x = 5;

	// time difference between two subsequent calls to the filter in seconds
	double dt;

	// dimension of radar measurement
	const unsigned nRadarDimension = 3;

	// dimension of laser measurement
	const unsigned nLaserDimension = 2;

	// declare measurement covariance matrix for radar
	MatrixXd S_Radar;

	// declare measurement covariance matrix for laser
	MatrixXd S_Laser;

	// declare a vector to hold radar measurement
	VectorXd actualRadarMeasurement;

	// declare a vector to hold laser measurement
	VectorXd actualLaserMeasurement;

	// declare a matrix to hold sigma points translated from state space to radar measurement space
	MatrixXd predictedSigmaPointsinRadarSpace;

	// declare a matrix to hold sigma points translated from state space to laser measurement space
	MatrixXd predictedSigmaPointsinLaserSpace;

	// declare a vector to hold estimated radar measurement based on predicted state
	VectorXd estimatedRadarMeasurement;

	// declare a vector to hold estimated laser measurement based on predicted state
	VectorXd estimatedLaserMeasurement;

	// declare a matrix to hold Kalman gain for Radar
	MatrixXd KalmanGainforRadar;

	// declare a matrix to hold Kalman gain for Laser
	MatrixXd KalmanGainforLaser;

	// instantiate tools
	Tools tools;

	/**
	 * Constructor
	 */
	UKF();

	/**
	 * Destructor
	 */
	virtual ~UKF();

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
};

#endif /* UKF_H */
