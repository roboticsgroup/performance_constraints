/*This is a demo of Performance Constraints for a robotic manipulator.

A description and the overall algorithm of the method is in:
- Dimeas, Fotios, Vassilis C. Moulianitis, and Nikos Aspragathos. 
"Manipulator performance constraints in human-robot cooperation." 
Robotics and Computer-Integrated Manufacturing (2017).

Author: Fotis Dimeas
Copyright 2022 Fotios Dimeas
 */

#include <cstdlib>
#include <iostream>

#include <armadillo> //Linear Algebra Library
#include <performanceConstraints.h>
#include <rclcpp/rclcpp.hpp>

using namespace std;

int main(int argc, char** argv) {
	cout << "Armadillo version: " << arma::arma_version::as_string() << endl;

	rclcpp::init(argc, argv);

	/*Initialize performance constraints
	* Arguments:
	* 1: w_cr for translation or combined
	* 2: w_th for translation or combined
	* 3: w_cr for rotation [OPTIONAL: Leave empty for combined indices]
	* 4: w_th for rotation [OPTIONAL: Leave empty for combined indices]
	* 5: lambda for translation or combined
	* 6: lambda for rotation  [OPTIONAL: Leave empty for combined indices]
	* 7: Performance indices [1, 2, 3]  1: Manipulability Index, 2: Minimum Singular Value, 3: Inverse Condition Number
	* 8: Calculation methods: _serial, _parallel
	* 9: Gradient with respect to: [_cartesian (default), _joints] Selectable for only gradient calculation
	*/

	/// Separate indices for translation or rotation
	PC pConstraints(	0.011,	0.03,	0.15,	0.5,	1.0,	1.0,	_manipulability,	_serial); //Using manipulability index
	// PC pConstraints(	0.03,	0.14,	0.1,	0.5,	1.0,	1.0,	_MSV,	_parallel); //Using MSV (better for human robot interaction)
	
	/// Combined indices
	// PC pConstraints(	0.01,	0.3,	1.0,	_MSV,	_serial); //Using MSV (better for human robot interaction)

	// Only calculate gradient
	// PC pConstraints(_manipulability,  _serial, _cartesian, 0, _JR_aug) ;

	pConstraints.setVerbose(0); //Set debug info. Comment or set to 0 to disable

	arma::vec q; 

	//for n=6 joints
	q = { 0., -M_PI/4., M_PI/2., 0., -M_PI/4., 0.0 }; //An example initial configuration 
	//for n=7 joints
	// q << 0. << -M_PI/4.  << 0 << M_PI/2. << 0. << -M_PI/4. << 0.0; //An example initial configuration 
	
	cout.precision(10); cout.setf(ios::fixed);
	
	cout << "Entering control loop...\n";
	arma::wall_clock timer; timer.tic();
	int i;
	for (i=0; i<100; i++) //this is supposed to be a simulation control loop
	{
		// Calculate PC in joint positions "q" [rad]
		pConstraints.updatePC(q); //Performance constraints are calculated in here
		arma::vec F = pConstraints.getSingularityTreatmentForce(); //

		// Calculate only the gradient without (no treatment force is calculated)
		// pConstraints.calculateGradient(q);
		// arma::vec A = pConstraints.getGradient();
		// A.t().print();
		// arma::vec A = pConstraints.getGradientScaled(10.0 * M_PI / 180.0, 5.0 * M_PI / 180.0);
		// F.t().print();

		if (pConstraints.checkForSingularity()) { //check for singularity since no handling of the robot's motion is done here
			cout << "Robot became singular. Stopping simulation after loops: " << i << endl; 
			break;	//stop the simulation
		}

		//gradually guide the robot to a singular configuration just to see the constraint forces
		q(1)+=0.01;
		q(3)-=0.01;
		q(5)+=0.01;
		
		//Put your controller here (e.g impedance or admittance)
		//...

		// cout << endl << endl;
	} 
	double time = timer.toc();
	cout << "Simulation completed in " <<  time<< "sec."<< endl;
	cout << "Average time per cycle: " << time/i << "sec."<< endl;

}