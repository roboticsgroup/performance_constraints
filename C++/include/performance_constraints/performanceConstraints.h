/*This is a demo of Performance Constraints for a robotic manipulator.

A description and the overall algorithm of the method is in:
- Dimeas, Fotios, Vassilis C. Moulianitis, and Nikos Aspragathos. 
"Manipulator performance constraints in human-robot cooperation." 
Robotics and Computer-Integrated Manufacturing (2017).

Author: Fotis Dimeas
Copyright 2022 Fotios Dimeas
 */


#ifndef PERFORMANCECONSTRAINTS_H
#define PERFORMANCECONSTRAINTS_H

#include "armadillo" //Linear Algebra Library
#include <thread>
#include <cstdlib>

#include <arl_core2/robot/robot_sim.hpp>
#include <arl_core2/robot/ros_model.hpp>
#include <rclcpp/rclcpp.hpp>

#ifndef M_PI 
#define M_PI 3.14159265359
#endif

/// Performance index
enum PerformanceIndex {
	_manipulability,
	_MSV,
	_iCN
};

/// Performance Constraints Calculation method
enum PCcalculation {
	_serial,
	_parallel
};

enum GradientWRT {
	_cartesian,
	_joints
};

enum JacobianToOptimize {
	_J,			//full jacobian
	_JT,		//position part only
	_JT_aug,	//augmented position part
	_JR,		//orientation part
	_JR_aug		//augmented orientaiton part
};

class PC
{
public:
	PC(double _crit_t, double _thres_t, double _crit_r, double _thres_r, double _lambda_t, double _lambda_r, PerformanceIndex _index, PCcalculation _method, bool _useAugmented = true);
	PC(double _crit_t, double _thres_t, double _lambda_t, PerformanceIndex _index, PCcalculation _method);
	PC(PerformanceIndex _index, PCcalculation _method, GradientWRT _gradient_type, bool _separate, JacobianToOptimize _optimJacobian = _J);
	~PC();
	void updatePC(const arma::vec q);
	void calculateGradient(const arma::vec q);
	double getSingularityTreatmentForce(int index){ return Favoid.at(index); } //return elements of the force/torque
	arma::vec getSingularityTreatmentForce() {return Favoid;}; //return the entire vector
	arma::vec getGradient() {return Aw;};
	arma::vec getGradientScaled(double clearance, double range);
	void setVerbose(int lvl) { verbose = lvl;}; //0 to disable
	int checkForSingularity();
	double getPerformanceIndex(int which=0);
	arma::vec getEquivalentStiffness();
	
private:	
	void init();
	void updatePC(); 
	void calculateGradient();
	void calcSingularityTreatmentForce(); //calculate the spring forces due to Singularity Treatment (call this from SingularityTreatment() )
	void findBestManip(PerformanceIndex option); //local search in the Cartesian tool frame (translations only) for best W (serial implementation)
	void findBestManip(); //Parallel performance constraints
	void threadpool_create(PerformanceIndex option); //Start thread-pool. Call this before entering control loop
	void threadpool_join(); //Joint threads on finish
	void Thread(int axis, PerformanceIndex option); //Thread implementation of Performance constraints
	
	void calcCurrentManipulability(); //calculate the manipulability index w=[T, R]
	double calcManipulability(const arma::mat J, const int T_R); //calculate and return the manipulability index
	double getManipulability(int T_R){ return w[T_R]; };  // Return the augmented translational or rotational manipulability index
	void calcCurrentSVD();
	double calcMSV(const arma::mat J, const int T_R);
	double calciCN(const arma::mat J, const int T_R);
	double getMSV(int T_R){ return msv[T_R]; }; // Return the current augmented translational or rotational Minimum Singular Value
	double getiCondNum(int T_R){ return icn[T_R]; }; // Return the current augmented translational or rotational inverse condition number
	
	arma::mat getJacobianToOptimize(arma::mat Jin);
	double calcCurrentIndex(arma::mat Jin);

	arma::mat get_Jsym_body(const arma::vec jntvalues); //Calculate the body Jacobian of the tool frame
	arma::mat get_Jsym_spatial(const arma::vec jntvalues); //Calculate spatial Jacobian of the base frame
	double getJointLimitScaling(double clearance, double range);

	void updateCurrentConfiguration(const arma::vec q); /** Update current robot configuration*/
	void chechInput(arma::vec q);

	arma::vec Aw; //index gradient
	int n; //number of joints used
	
	//Performance constraints
	double K_w; //Singularity treatment spring 
	double K_rot; //Singularity treatment spring (for rotations)
	arma::vec::fixed<2> ST_current; //Singularity treatment current value (e.g. manipulability index)
	arma::vec::fixed<6> Favoid, dp;
	double dx, grad_w;
	arma::vec Qv, Qinit, Q_measured, dq, qlim; 
	std::vector<std::thread> pconstraints;
	int stop_pConstraints_pool;
	arma::vec updateConstraints;
	double crit_t, thres_t, crit_r, thres_r, lambda_t, lambda_r;
	PerformanceIndex PC_index; //performance constraints option
	PCcalculation PC_Calc_Method; //serial or parallel calculation
	GradientWRT gradient_type; //with respect to cartesian frame or joints
	JacobianToOptimize optimJacobian;

	arma::mat J; //This is the Jacobian matrix at the current configuration of the robot. It is calculated based on the loaded URDF using q.
	arma::mat J_sym; //The Jacobian matrix in neighbor configurations based on the loaded URDF 
	
	// 2 part indices [translational, rotational]
	arma::vec::fixed<2> w; //manipulability index of current position, calculated in every loop inside msrJacobian()
	arma::vec::fixed<2> msv; //Minimum singular value
	arma::vec::fixed<2> icn; //inverse condition number
	
	int verbose;
	arma::wall_clock timer;
	bool separate, useAugmented;
	double qlim_max;

	std::shared_ptr<arl::robot::Robot> sim_robot;
	std::shared_ptr<arl::robot::ROSModel> model;
};

#endif 