/*This is a demo of Performance Constraints for a robotic manipulator.
The KUKA LWR 4+ is provided as an example here.

A description and the overall algorithm of the method is in:
- Dimeas, Fotios, Vassilis C. Moulianitis, and Nikos Aspragathos. 
"Manipulator performance constraints in human-robot cooperation." 
Robotics and Computer-Integrated Manufacturing (2017).

Author: Fotis Dimeas
Copyright 2020 Fotios Dimeas
 */

 #include "performanceConstraints.h"

/*! \brief Initialize the Performance Constraints
*
* Call the constructor before entering the control loop for separate indices.
* \param _method Available calculation methods: _serial, _parallel, _parallel_nonblock
*/
PC::PC(	double _crit_t, double _thres_t, 
	double _crit_r, double _thres_r, 
	double _lambda_t, double _lambda_r, 
	int _index, 
	PCcalculation _method) {

	crit_t 		= _crit_t;
	thres_t 	= _thres_t;
	crit_r 		= _crit_r;
	thres_r 	= _thres_r;
	lambda_t 	= _lambda_t;
	lambda_r 	= _lambda_r;
	PC_index 	= _index;
	PC_Calc_Method	= _method;

	separate 	= true;

	init();
}

/* Call the constructor before entering the control loop for combined indices.
* \param _method Available calculation methods: _serial, _parallel, _parallel_nonblock
*/
PC::PC(	double _crit_t, double _thres_t, 
	double _lambda_t,  
	int _index, 
	PCcalculation _method) {

	crit_t 		= _crit_t;
	thres_t 	= _thres_t;
	lambda_t 	= _lambda_t;
	PC_index 	= _index;
	PC_Calc_Method	= _method;

	separate 	= false;
	
	init();
}

void PC::init() {
	if (PC_index<1 || PC_index>3) {
		std::cout << "Unknown calculation option '" << PC_index << "'. Use between 1 for manipulability, 2 for MSV or 3 for iCN." << std::endl;
		std::cout << "Using 2 as a defult" << std::endl;
		PC_index = 2;
	}

	n=7; //by default use all 7 joints for differential inverse kinematics etc.
	K_w = 0.0;
	K_rot = 0.0;
	Aw.fill(1.0);
	Favoid.fill(0.0);
	dx = 1e-5; //infinitesimal movement
	w_all_prev = 0.0;

	if (PC_Calc_Method != _serial) {
		threadpool_create(PC_index); //start thread pool
	}

	
	std::cout << "PC has been initialized using method: "; 
	switch (PC_index){
		case 1:			//Manipulability index
			std::cout << "'Manipulability Index' ";
			break;
		case 2:			//Minimum singulr value
			std::cout << "'Minimum Singular Value' ";
			break;
		case 3:			//Condition number
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
		w(0) = sqrt(det( J.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(3,5))*J.rows(3,5)) * J.rows(0,2).t() )); //translational strong sense
		w(1) = sqrt(det( J.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(0,2))*J.rows(0,2)) * J.rows(3,5).t() )); //rotational strong sense	
	}
	else {
		w(0) = sqrt(det( J * J.t() )); 
		w(1) = w(0); //only used in the case of seperate=false
	}
}

/// Calculate and return the augmented transl. or rotational manipulability
double PC::calcManipulability(const arma::mat J, const int T_R)
{
	if (separate) {
		if (T_R==0) { //translational
			// return sqrt(det(J.rows(0,2)*J.rows(0,2).t())); //weak sense
			return sqrt(det( J.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(3,5))*J.rows(3,5)) * J.rows(0,2).t() )); //strong sense
		}
		else { //rotational
			// return sqrt(det(J.rows(3,5)*J.rows(3,5).t())); //weak sense
			return sqrt(det( J.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(0,2))*J.rows(0,2)) * J.rows(3,5).t() )); //strong sense
		}
	}
	else {
		return sqrt(det( J * J.t() )); 
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
double PC::calcMSV(const arma::mat J, const int T_R) {
	arma::vec sigma;
	if (separate) {
		if (T_R==0) { //translational
			sigma = arma::svd( J.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(3,5))*J.rows(3,5)) ); 
		}
		else { //rotational
			sigma = arma::svd( J.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(0,2))*J.rows(0,2)) ); 
		}
	}
	else {
		sigma = arma::svd( J );
	}
	return arma::min(sigma);
}

/// Calculate and return the augmented transl. or rotational inverse condition number (LCI)
double PC::calciCN(const arma::mat J, const int T_R) {
	arma::vec sigma;
	if (separate) {
		if (T_R==0) { //translational
			sigma = arma::svd( J.rows(0,2) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(3,5))*J.rows(3,5)) ); 
		}
		else { //rotational
			sigma = arma::svd( J.rows(3,5) * (arma::eye<arma::mat>(n,n)-arma::pinv(J.rows(0,2))*J.rows(0,2)) ); 
		}
	}
	else {
		sigma = arma::svd( J );
	}
	return arma::min(sigma)/arma::max(sigma);
}

/*! \brief Calculate the Performance Constraints
*
* Call this function from the control loop.
* If the current performance value drops below the threshold, the grad is calculated.
*/
void PC::updatePC()
{
	if (verbose) timer.tic();
	
	if (PC_index==1)
		calcCurrentManipulability(); 
	else
		calcCurrentSVD();  

	//Select performance index
	switch (PC_index){
		case 1:			//Manipulability index
			ST_current = w;
			break;
		case 2:			//Minimum singulr value
			ST_current = msv;
			break;
		case 3:			//Condition number
			ST_current = icn;
			break;
	}

	//Select calculation method (Serial or parallel)
	if (PC_Calc_Method == _serial)
		findBestManip(PC_index); //serial calculation
	else
		findBestManip(); //parallel calculation

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
		if (separate) std::cout << "Current performance (Lin,Rot): " << ST_current(0) << " | " << ST_current(1) << std::endl;
		else std::cout << "Current performance (combined): " << ST_current(0) << std::endl;
		std::cout << "Constraint forces: "; Favoid.t().raw_print();
		std::cout << "Took: " << time << "sec" << std::endl;
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
* \param option selects between 1:manipulability, 2:minimum singular value, 3:Inv condition number
*/
void PC::findBestManip(int option){
	Qinit = Q_measured; //copy current measured joint values (q_0)
	
	for (int axis=0; axis<6; axis++){ //for each Cartesian direction
		int T_R = (axis<3) ? 0 : 1; //0 if translational, 1 if rotational [If seperate=false, the value of T_R doesn't play any role]
		
		dp.fill(0.0); 
		dp.at(axis)=dx; //virtual Cart. velocity for local search
		
		Qv = arma::pinv(J)*dp + Qinit; //new virtual joint values
		
		get_Jsym(Qv, J_sym); //calculate J for those new virtual joint values (This overwrites the J7_sym global variable)
		
		switch (option){
			case 1:			//Manipulability metric
				grad_w = calcManipulability(J_sym, T_R)-w[T_R];
				break;
			case 2:			//Minimum singulr value metric
				grad_w = calcMSV(J_sym, T_R) - msv[T_R];
				break;
			case 3:			//Minimum singulr value metric
				grad_w = calciCN(J_sym, T_R) - icn[T_R];
				break;
		}

		Aw.at(axis) = grad_w / dx; 
	}
}

/*! \brief Calculate the gradient of the performance index wrt to the Cartesian frame. Parallel (non-blocking) implementation
*
* Non-blocking means that in each loop the constraints from the previous iteration will be applied
* The option selects between 1:manipulability, 2:minimum singular value, 3:Inv condition number
*/
void PC::findBestManip(){
	Qinit = Q_measured; //copy current measured joint values (q_0)
	
	for (int axis=0; axis<6; axis++) 
		updateConstraints[axis] = 1; //send signal at the thread pool to update measurements
	
	//join (all of them must become 0 to continue) [becomes non-blocking otherwise]
	if (PC_Calc_Method == _parallel) {
		while (	updateConstraints[0]==1 || 
				updateConstraints[1]==1 || 
				updateConstraints[2]==1 || 
				updateConstraints[3]==1 || 
				updateConstraints[4]==1 || 
				updateConstraints[5]==1 	) {
			usleep(5); //sleep for a while to avoid CPU 100%
		}
	}
}

/*! \brief Thread impelemtation of performance constraints
* 
* Thread-pool like implemetation. Start a thread for each of the Cartesian directions with the desired option (index).
* More computationally heavy because of matrix copies and multiple declarations but non-blocking
*/
void PC::Thread(int axis, int option){
	int T_R = (axis<3) ? 0 : 1; //0 if translational, 1 if rotational [If seperate=false, the value of T_R doesn't play any role]
	double grad_w;
	arma::mat Jc(6,n);

	while(!stop_pConstraints_pool) {
		if (updateConstraints[axis]) {

			arma::vec dp(6); dp.fill(0.0); 
			dp.at(axis)=dx; //virtual Cart. velocity for local search
			
			arma::vec Qv = arma::pinv(J)*dp + Qinit; //new virtual joint values
			
			get_Jsym(Qv, Jc); //calculate J for those new virtual joint values 
			
			switch (option){
				case 1:			//Manipulability metric
					grad_w = calcManipulability(Jc, T_R)-w[T_R];
					break;
				case 2:			//Minimum singulr value metric
					grad_w = calcMSV(Jc, T_R) - msv[T_R];
					break;
				case 3:			//Minimum singulr value metric
					grad_w = calciCN(Jc, T_R) - icn[T_R];
					break;
			}
			
			Aw.at(axis) = grad_w / dx; 

			updateConstraints[axis] = 0; //finished
		}
		else
			usleep(10); //sleep for a while to avoid CPU 100%
	}
}

/*! \brief Spawn threads for parallel calculation of performance constraints
*
*/
void PC::threadpool_create(int option) {
	if (verbose) std::cout << "Creating thread pool...";
	stop_pConstraints_pool = 0;
	for (int axis=0; axis<6; axis++)
		updateConstraints[axis] = 0;

	for (int axis=0; axis<6; axis++){ //for each Cartesian direction
		pconstraints[axis] = std::thread(&PC::Thread, this, axis, option);
	}
	if (verbose) std::cout << " done!" << std::endl;
}

/*! \brief Joint threads of parallel calculation of performance constraints
*
*/
void PC::threadpool_join() {
	stop_pConstraints_pool = 1;
	for (int axis=0; axis<6; axis++){ //for each Cartesian direction
		pconstraints[axis].join();
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