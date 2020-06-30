/*This is a demo of Performance Constraints for a robotic manipulator.
The KUKA LWR 4+ is provided as an example here.

A description and the overall algorithm of the method is in:
- Dimeas, Fotios, Vassilis C. Moulianitis, and Nikos Aspragathos. 
"Manipulator performance constraints in human-robot cooperation." 
Robotics and Computer-Integrated Manufacturing (2017).

Author: Fotis Dimeas
Copyright 2020 Fotios Dimeas
 */

#include <performance_constraints/performanceConstraints.h>

/*! \brief Initialize the Performance Constraints
*
* Call the constructor before entering the control loop for separate indices.
* \param _method Available calculation methods: _serial, _parallel, _parallel_nonblock
*/
PC::PC(	double _crit_t, double _thres_t, 
	double _crit_r, double _thres_r, 
	double _lambda_t, double _lambda_r, 
	PerformanceIndex _index, 
	PCcalculation _method) {

	crit_t 		= _crit_t;
	thres_t 	= _thres_t;
	crit_r 		= _crit_r;
	thres_r 	= _thres_r;
	lambda_t 	= _lambda_t;
	lambda_r 	= _lambda_r;
	PC_index 	= _index;
	PC_Calc_Method	= _method;

	gradient_type = _cartesian;
	separate 	= true; //different thresholds have been provided for each

	init();
}

/* Call the constructor before entering the control loop for combined indices in position and orientation.
* \param _method Available calculation methods: _serial, _parallel, _parallel_nonblock
*/
PC::PC(	double _crit_t, double _thres_t, 
	double _lambda_t,  
	PerformanceIndex _index, 
	PCcalculation _method) {

	crit_t 		= _crit_t;
	thres_t 	= _thres_t;
	lambda_t 	= _lambda_t;
	PC_index 	= _index;
	PC_Calc_Method	= _method;;

	gradient_type = _cartesian;
	separate 	= false; //only one set of thresholds has been provided for the combined index
	
	init();
}

/* Call this constructor to calculate only the Gradient.
* \param _method Available calculation methods: _serial, _parallel, _parallel_nonblock
*/
PC::PC(	PerformanceIndex _index, 
		PCcalculation _method,
		GradientWRT _gradient_type,
		bool _separate) {

	PC_index 	= _index;
	PC_Calc_Method	= _method;
	gradient_type = _gradient_type;

	if (gradient_type == _joints) {
		if (verbose) std::cout << "Warning! Using combined indices when the gradient is wrt. the joints." << std::endl;
		separate = false; //when grad is wrt. joints, there is no point to separate
	}
	else if (gradient_type == _cartesian)
		separate = _separate; 
	else {
		std::cout << "Error! Wrong gradient calculation type" << std::endl;
		return;
	}
	
	init();
}

void PC::init() {
	if ( !(PC_index==_manipulability || PC_index==_MSV || PC_index==_iCN)) {
		std::cout << "Unknown calculation option '" << PC_index << "'. Use: _manipulability, _MSV or _iCN." << std::endl;
		std::cout << "Using _MSV as a defult" << std::endl;
		PC_index = _MSV;
	}

	n=7; //by default use all 7 joints for differential inverse kinematics etc.
	
	K_w = 0.0;
	K_rot = 0.0;
	
	if (gradient_type == _cartesian) {
		Aw = arma::ones<arma::vec>(6);
		updateConstraints =  arma::zeros<arma::vec>(6);
	}
	else { 
		Aw = arma::ones<arma::vec>(n);
		updateConstraints =  arma::zeros<arma::vec>(n);
	}

	dq = arma::ones<arma::vec>(n);
	Favoid.fill(0.0);

	dx = 1e-5; //infinitesimal movement

	if (PC_Calc_Method != _serial) {
		threadpool_create(PC_index); //start thread pool
	}

	std::cout << "PC has been initialized using method: "; 
	switch (PC_index){
		case _manipulability:			//Manipulability index
			std::cout << "'Manipulability Index' ";
			break;
		case _MSV:			//Minimum singulr value
			std::cout << "'Minimum Singular Value' ";
			break;
		case _iCN:			//Inverse condition number
			std::cout << "'Inverse Condition Number' ";
			break;
	}
	std::cout << std::endl;

	verbose = 0; //default
}

/*! \brief Destructor
*/
PC::~PC() {
	if (PC_Calc_Method!=_serial) { //only in parallel mode
		threadpool_join();
		if (verbose) std::cout << "Threads have been shut down safely." << std::endl;
	}
	else
		if (verbose) std::cout << "Destroying PC object..." << std::endl;
	
	if (verbose) std::cout << "All done!" << std::endl;
}

/// Calculate the current augmented translational and rotational manipulability from J (Yoshikawa, 1990)
void PC::calcCurrentManipulability() {
	if (separate) {
		w(0) = sqrt(arma::det( J.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(3,5))*J.rows(3,5)) * J.rows(0,2).t() )); //translational strong sense
		w(1) = sqrt(arma::det( J.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(0,2))*J.rows(0,2)) * J.rows(3,5).t() )); //rotational strong sense	
	}
	else {
		w(0) = sqrt(arma::det( J * J.t() )); 
		w(1) = w(0); //only used in the case of seperate=false
	}
	// if (verbose) std::cout << "Calculated current manipulability: " << w(0) << " | " << w(1) << std::endl;
}

/// Calculate and return the augmented transl. or rotational manipulability
double PC::calcManipulability(const arma::mat J_, const int T_R)
{
	if (separate) {
		if (T_R==0) { //translational
			// return sqrt(arma::det(J_.rows(0,2)*J_.rows(0,2).t())); //weak sense
			return sqrt(arma::det( J_.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J_.rows(3,5))*J_.rows(3,5)) * J_.rows(0,2).t() )); //strong sense
		}
		else { //rotational
			// return sqrt(arma::det(J_.rows(3,5)*J_.rows(3,5).t())); //weak sense
			return sqrt(arma::det( J_.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J_.rows(0,2))*J_.rows(0,2)) * J_.rows(3,5).t() )); //strong sense
		}
	}
	else {
		return sqrt(arma::det( J_ * J_.t() )); 
	}
}

/** Calculate SVD on the current Jacobian and return MSV and CN.*/
void PC::calcCurrentSVD(){
	if (separate) {
		//translational
		arma::vec sigmaT = arma::svd( J.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(3,5))*J.rows(3,5)) ); 
		msv(0) = arma::min(sigmaT);
		icn(0) = arma::min(sigmaT)/arma::max(sigmaT);

		//rotational
		arma::vec sigmaR = arma::svd( J.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(0,2))*J.rows(0,2)) ); 
		msv(1) = arma::min(sigmaR);
		icn(1) = arma::min(sigmaR)/arma::max(sigmaR);
	}
	else {
		arma::vec sigma = arma::svd( J );
		msv(0) = arma::min(sigma);
		icn(0) = arma::min(sigma)/arma::max(sigma);
		msv(1) = msv(0); //only used in the case of seperate=false
		icn(1) = icn(0); //only used in the case of seperate=false
	}
}

/// Calculate and return the augmented transl. or rotational MSV
double PC::calcMSV(const arma::mat J_, const int T_R) {
	arma::vec sigma;
	if (separate) {
		if (T_R==0) { //translational
			sigma = arma::svd( J_.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J_.rows(3,5))*J_.rows(3,5)) ); 
		}
		else { //rotational
			sigma = arma::svd( J_.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J_.rows(0,2))*J_.rows(0,2)) ); 
		}
	}
	else {
		sigma = arma::svd( J_ );
	}
	return arma::min(sigma);
}

/// Calculate and return the augmented transl. or rotational inverse condition number (LCI)
double PC::calciCN(const arma::mat J_, const int T_R) {
	arma::vec sigma;
	if (separate) {
		if (T_R==0) { //translational
			sigma = arma::svd( J_.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J_.rows(3,5))*J_.rows(3,5)) ); 
		}
		else { //rotational
			sigma = arma::svd( J_.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J_.rows(0,2))*J_.rows(0,2)) ); 
		}
	}
	else {
		sigma = arma::svd( J_ );
	}
	return arma::min(sigma)/arma::max(sigma);
}

/*! \brief Calculate the Performance Constraints
*
* Call this function from outside, simply by providing the current q.
*/
void PC::updatePC(const arma::vec q){
	updateCurrentConfiguration(q); //measure the robot's configuration and put it here
	
	J = get_Jsym_spatial(q); //calculate J from current q

	updatePC();

}

/*! \brief Calculate the Performance Constraints
*
* Call this function from outside, but fist need to updateCurrentConfiguration(q), get_Jsym_spatial(q, J) and updateCurrentJacobian(J)
* If the current performance value drops below the threshold, the grad is calculated.
*/
void PC::updatePC()
{
	if (verbose) timer.tic();
	
	calculateGradient(); //The point of PC is to work with the _cartesian gradient 

	//Calculate metrics
	if (ST_current(0) < thres_t) { //translation	
		K_w = -lambda_t*( 1.0/(ST_current(0)-crit_t)-1.0/(thres_t-crit_t) ); //asymptotical increase
	}
	else {
		K_w = 0.0; //spring force disabled
	}

	if (separate) {
		if (ST_current(1) < thres_r) { //rotation
			K_rot = -lambda_r*( 1.0/(ST_current(1)-crit_r)-1.0/(thres_r-crit_r) ); //asymptotical increase
		}
		else {
			K_rot = 0.0; //spring force disabled
		}
	}

	//Calculate the forces/torques
	calcSingularityTreatmentForce(); 

	if (verbose) {
		double time = timer.toc();
		std::cout << "Constraint forces: "; Favoid.t().raw_print();
		std::cout << "Took: " << time << "sec" << std::endl;
	}
}

/// Wrapper of calculateGradient() that gets the current joint positions as an argument.
void PC::calculateGradient(const arma::vec q){
	updateCurrentConfiguration(q); //measure the robot's configuration and put it here
	
	J = get_Jsym_spatial(q); //calculate J from current q

	calculateGradient();

}

///  To call this function from outside, make sure you call first: updateCurrentConfiguration(q), get_Jsym_spatial(q, J) and updateCurrentJacobian(J)
void PC::calculateGradient() {
	//Select performance index
	switch (PC_index){
		case _manipulability:			//Manipulability index
			calcCurrentManipulability(); 
			ST_current = w;
			break;
		case _MSV:			//Minimum singulr value
			calcCurrentSVD(); 
			ST_current = msv;
			break;
		case _iCN:			//Condition number
			calcCurrentSVD(); 
			ST_current = icn;
			break;
	}

	//Select calculation method (Serial or parallel)
	if (PC_Calc_Method == _serial)
		findBestManip(PC_index); //serial calculation
	else
		findBestManip(); //parallel calculation

	if (verbose) {
		if (separate) std::cout << "Current performance (Lin,Rot): " << ST_current(0) << " | " << ST_current(1) << std::endl;
		else std::cout << "Current performance (combined): " << ST_current(0) << std::endl;
	}
}

/*! \brief Calculate the reaction force from the Performance constraints
*
* Calculate the reaction force in a certain Cartesian axis.
* The reaction force is the product of the distance from the performance threshold and 
* the grad of the performance index with respect to the Cartesian axis.
*/
void PC::calcSingularityTreatmentForce(){
	for (int index=0; index<6; index++) { //for each axis
		if (separate) {
			if (index<=2) { //for translational directions
				Favoid.at(index) = - Aw.at(index)*K_w;
			}
			else if	(index<=5) { //for rotational directions
				Favoid.at(index) = - Aw.at(index)*K_rot;
			}
		}
		else { //for combined index
			Favoid.at(index) = - Aw.at(index)*K_w;
		}
		
	}
}

/*! \brief Calculate the gradient of the performance index wrt to the Cartesian frame. Serial Implementation
*
* \param option selects between _manipulability, _MSV, _iCN
*/
void PC::findBestManip(PerformanceIndex option){
	Qinit = Q_measured; //copy current measured joint values (q_0)
	
	for (int axis=0; axis<Aw.size(); axis++){ //for each direction (cartesian or joint)
		int T_R;
		if (gradient_type == _cartesian) {
			T_R = (axis<3) ? 0 : 1; //0 if translational, 1 if rotational [If seperate=false, the value of T_R doesn't play any role]
			dp.fill(0.0); 
			dp.at(axis)=dx; //virtual Cart. velocity for local search
			
			Qv = arma::pinv(J)*dp + Qinit; //new virtual joint values
		}
		else if (gradient_type == _joints) {
			T_R = 0;
			dq.fill(0.0); 
			dq.at(axis)=dx; //infinitesimal joint movement
			
			Qv = dq + Qinit; //new virtual joint values
		}
		
		J_sym = get_Jsym_spatial(Qv); //calculate J for those new virtual joint values (This overwrites the J_sym global variable)
		
		switch (option){
			case _manipulability:			//Manipulability metric
				grad_w = calcManipulability(J_sym, T_R)-w[T_R];
				break;
			case _MSV:			//Minimum singulr value metric
				grad_w = calcMSV(J_sym, T_R) - msv[T_R];
				break;
			case _iCN:			//Minimum singulr value metric
				grad_w = calciCN(J_sym, T_R) - icn[T_R];
				break;
		}

		Aw.at(axis) = grad_w / dx; 
	}	
}

/*! \brief Calculate the gradient of the performance index wrt to the Cartesian frame. Parallel implementation
*
* _parallel is a blocking function
* _parallel_nonblock means that in each loop, the last calculated constraints from the previous iteration will be applied [WARNING! Not Advised!]
*/
void PC::findBestManip(){
	Qinit = Q_measured; //copy current measured joint values (q_0)
	
	// for (int axis=0; axis<6; axis++) 
	// 	updateConstraints[axis] = 1; //send signal at the thread pool to update measurements
	updateConstraints.ones();

	//join (all of them must become 0 to continue) [becomes non-blocking otherwise]
	if (PC_Calc_Method == _parallel) {
		while (	any(updateConstraints) ) {
			// usleep(1); //sleep for a while to avoid CPU 100%
		}
	}

}

/*! \brief Thread impelemtation of performance constraints
* 
* Thread-pool like implemetation. Start a thread for each of the Cartesian directions with the desired option (index).
* More computationally heavy because of matrix copies and multiple declarations but non-blocking
*/
void PC::Thread(int axis, PerformanceIndex option){
	int T_R;
	double grad_w;
	arma::vec Qv(n), dp, dq;
	if (gradient_type == _cartesian) {
		T_R = (axis<3) ? 0 : 1; //0 if translational, 1 if rotational [If seperate=false, the value of T_R doesn't play any role]
		dp.resize(6);
	}
	else if (gradient_type == _joints) {
		T_R = 0;
		dq.resize(n);
	}
	arma::mat Jc(6,n);

	while(!stop_pConstraints_pool) {
		if (updateConstraints(axis)) {

			if (gradient_type == _cartesian) {
				dp.fill(0.0); 
				dp.at(axis)=dx; //virtual Cart. velocity for local search
				
				Qv = arma::pinv(J)*dp + Qinit; //new virtual joint values
			}
			else if (gradient_type == _joints) {
				// arma::vec dq(n);
				dq.fill(0.0); 
				dq.at(axis)=dx; //infinitesimal joint movement
				
				Qv = dq + Qinit; //new virtual joint values
			}
			
			Jc = get_Jsym_spatial(Qv); //calculate J for those new virtual joint values 
			
			switch (option){
				case _manipulability:			//Manipulability metric
					grad_w = calcManipulability(Jc, T_R)-w[T_R];
					break;
				case _MSV:			//Minimum singulr value metric
					grad_w = calcMSV(Jc, T_R) - msv[T_R];
					break;
				case _iCN:			//Minimum singulr value metric
					grad_w = calciCN(Jc, T_R) - icn[T_R];
					break;
			}
			
			Aw.at(axis) = grad_w / dx; 

			updateConstraints(axis) = 0.; //finished
		}
		else
			usleep(1); //sleep for a while to avoid CPU 100%. Sleep seems to be necessary within the thread.
	}
}

/*! \brief Spawn threads for parallel calculation of performance constraints
*
*/
void PC::threadpool_create(PerformanceIndex option) {
	if (verbose) std::cout << "Creating thread pool...";
	stop_pConstraints_pool = 0;

	updateConstraints.zeros();

	for (int axis=0; axis<updateConstraints.size(); axis++){ //for each Cartesian direction
		pconstraints.push_back(std::thread(&PC::Thread, this, axis, option));
	}
	if (verbose) std::cout << " done! " << pconstraints.size() << " threads created." << std::endl;
}

/*! \brief Joint threads of parallel calculation of performance constraints
*
*/
void PC::threadpool_join() {
	stop_pConstraints_pool = 1;
	for (int axis=0; axis<updateConstraints.size(); axis++){ //for each Cartesian direction
		pconstraints.at(axis).join();
	}
}

/*! \brief Update current robot configuration
*
*/
void PC::updateCurrentConfiguration(const arma::vec q){
	Q_measured = q;
}

/*! \brief Update current robot Jacobian
*
*/
void PC::updateCurrentJacobian(const arma::mat J_current){
	J = J_current;
}

double PC::getPerformanceIndex(int which){
	return ST_current(which);
}

int PC::checkForSingularity() {
	if (separate) {
		if (ST_current(0) < crit_t || ST_current(1) < crit_r) 
			return 1; //means that you should be carefull
		else 
			return 0;
	}
	else {
		if (ST_current(0) < crit_t) 
			return 1; //means that you should be carefull
		else 
			return 0;
	}
}