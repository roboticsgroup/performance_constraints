/*This is a demo of Performance Constraints for a robotic manipulator.
The KUKA LWR 4+ is provided as an example here.

A description and the overall algorithm of the method is in:
- Dimeas, Fotios, Vassilis C. Moulianitis, and Nikos Aspragathos. 
"Manipulator performance constraints in human-robot cooperation." 
Robotics and Computer-Integrated Manufacturing (2017).

In here we calculate the Jacobian matrix (6x7) for the KUKA LWR manipulator (expressed at the tool frame)
The calculations are derived from Matlab symbolic toolbox and Ccode creator.

Author: Fotis Dimeas
Copyright 2017 Fotios Dimeas
*/

#include "performanceConstraints.h"

void PC::get_Jsym(const arma::vec jntvalues, arma::mat& J) {
	
	double t1 = jntvalues.at(0);
	double t2 = jntvalues.at(1);
	double t3 = jntvalues.at(2);
	double t4 = jntvalues.at(3);
	double t5 = jntvalues.at(4);
	double t6 = jntvalues.at(5);
	double t7 = jntvalues.at(6);

	J.at(0,0) = (sin(t1)*(sin(t2)*(sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))+cos(t4)*cos(t7)*sin(t6))+cos(t2)*(cos(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+sin(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))))+cos(t1)*(sin(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))-cos(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))))*(cos(t1)*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0))-sin(t1)*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2)))-(cos(t1)*(sin(t2)*(sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))+cos(t4)*cos(t7)*sin(t6))+cos(t2)*(cos(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+sin(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))))-sin(t1)*(sin(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))-cos(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))))*(cos(t1)*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t1)*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0)));
	J.at(0,1) = -(sin(t2)*(sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))+cos(t4)*cos(t7)*sin(t6))+cos(t2)*(cos(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+sin(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))))*(sin(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+cos(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0))-(cos(t2)*(sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))+cos(t4)*cos(t7)*sin(t6))-sin(t2)*(cos(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+sin(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))))*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0));
	J.at(0,2) = -(cos(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+sin(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5)))*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+(sin(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))-cos(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5)))*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2));
	J.at(0,3) = (sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))+cos(t4)*cos(t7)*sin(t6))*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2));
	J.at(0,4) = cos(t5)*sin(t6)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5))*(3.9E1/5.0E2)+sin(t5)*sin(t6)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))*(3.9E1/5.0E2);
	J.at(0,5) = pow(cos(t6),2.0)*cos(t7)*(3.9E1/5.0E2)+cos(t7)*pow(sin(t6),2.0)*(3.9E1/5.0E2);
	J.at(1,0) = -(cos(t1)*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t1)*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0)))*(cos(t1)*(sin(t2)*(sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))-cos(t4)*sin(t6)*sin(t7))+cos(t2)*(cos(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+sin(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7))))-sin(t1)*(sin(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))-cos(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7))))+(cos(t1)*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0))-sin(t1)*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2)))*(cos(t1)*(sin(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))-cos(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7)))+sin(t1)*(sin(t2)*(sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))-cos(t4)*sin(t6)*sin(t7))+cos(t2)*(cos(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+sin(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7)))));
	J.at(1,1) = -(sin(t2)*(sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))-cos(t4)*sin(t6)*sin(t7))+cos(t2)*(cos(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+sin(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7))))*(sin(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+cos(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0))-(cos(t2)*(sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))-cos(t4)*sin(t6)*sin(t7))-sin(t2)*(cos(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+sin(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7))))*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0));
	J.at(1,2) = -(cos(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+sin(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7)))*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+(sin(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))-cos(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7)))*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2));
	J.at(1,3) = (cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2))+(sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))-cos(t4)*sin(t6)*sin(t7))*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2));
	J.at(1,4) = cos(t5)*sin(t6)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7))*(3.9E1/5.0E2)+sin(t5)*sin(t6)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))*(3.9E1/5.0E2);
	J.at(1,5) = pow(cos(t6),2.0)*sin(t7)*(-3.9E1/5.0E2)-pow(sin(t6),2.0)*sin(t7)*(3.9E1/5.0E2);
	J.at(2,0) = (cos(t1)*(sin(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))-cos(t3)*sin(t5)*sin(t6))+sin(t1)*(cos(t2)*(cos(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+sin(t3)*sin(t5)*sin(t6))-sin(t2)*(cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6))))*(cos(t1)*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0))-sin(t1)*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2)))+(sin(t1)*(sin(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))-cos(t3)*sin(t5)*sin(t6))-cos(t1)*(cos(t2)*(cos(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+sin(t3)*sin(t5)*sin(t6))-sin(t2)*(cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6))))*(cos(t1)*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t1)*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0)));
	J.at(2,1) = (sin(t2)*(cos(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+sin(t3)*sin(t5)*sin(t6))+cos(t2)*(cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6)))*(cos(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))-sin(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0))-(cos(t2)*(cos(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+sin(t3)*sin(t5)*sin(t6))-sin(t2)*(cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6)))*(sin(t2)*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+cos(t2)*(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2)+2.0/5.0));
	J.at(2,2) = -(cos(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+sin(t3)*sin(t5)*sin(t6))*(sin(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))-cos(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2))+(sin(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))-cos(t3)*sin(t5)*sin(t6))*(cos(t3)*(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))+sin(t3)*sin(t5)*sin(t6)*(3.9E1/5.0E2));
	J.at(2,3) = -(sin(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)-cos(t4)*cos(t5)*sin(t6)*(3.9E1/5.0E2))*(cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6))+(cos(t4)*(cos(t6)*(3.9E1/5.0E2)+2.0/5.0)+cos(t5)*sin(t4)*sin(t6)*(3.9E1/5.0E2))*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6));
	J.at(3,0) = -cos(t2)*(sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))+cos(t4)*cos(t7)*sin(t6))+sin(t2)*(cos(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+sin(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5)));
	J.at(3,1) = -sin(t3)*(cos(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t7)*sin(t4)*sin(t6))+cos(t3)*(cos(t5)*sin(t7)+cos(t6)*cos(t7)*sin(t5));
	J.at(3,2) = -sin(t4)*(sin(t5)*sin(t7)-cos(t5)*cos(t6)*cos(t7))-cos(t4)*cos(t7)*sin(t6);
	J.at(3,3) = -cos(t5)*sin(t7)-cos(t6)*cos(t7)*sin(t5);
	J.at(3,4) = -cos(t7)*sin(t6);
	J.at(3,5) = sin(t7);
	J.at(4,0) = -cos(t2)*(sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))-cos(t4)*sin(t6)*sin(t7))+sin(t2)*(cos(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+sin(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7)));
	J.at(4,1) = -sin(t3)*(cos(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+sin(t4)*sin(t6)*sin(t7))+cos(t3)*(cos(t5)*cos(t7)-cos(t6)*sin(t5)*sin(t7));
	J.at(4,2) = -sin(t4)*(cos(t7)*sin(t5)+cos(t5)*cos(t6)*sin(t7))+cos(t4)*sin(t6)*sin(t7);
	J.at(4,3) = -cos(t5)*cos(t7)+cos(t6)*sin(t5)*sin(t7);
	J.at(4,4) = sin(t6)*sin(t7);
	J.at(4,5) = cos(t7);
	J.at(5,0) = sin(t2)*(cos(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+sin(t3)*sin(t5)*sin(t6))+cos(t2)*(cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6));
	J.at(5,1) = -sin(t3)*(cos(t6)*sin(t4)-cos(t4)*cos(t5)*sin(t6))+cos(t3)*sin(t5)*sin(t6);
	J.at(5,2) = cos(t4)*cos(t6)+cos(t5)*sin(t4)*sin(t6);
	J.at(5,3) = -sin(t5)*sin(t6);
	J.at(5,4) = cos(t6);
	J.at(5,6) = 1.0;

	}