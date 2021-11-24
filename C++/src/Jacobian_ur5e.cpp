/*Symbolic Jacobian for UR5e 

Author: Fotis Dimeas
Copyright 2021 Fotios Dimeas
*/

#include <performance_constraints/performanceConstraints.h>

//Jacobian expressed on the tool frame
arma::mat PC::get_Jsym_body(const arma::vec jntvalues) {
	
	arma::mat J(6,6);
	J.zeros(); //very important to initialize to zero!

	double t1 = jntvalues.at(0);
	double t2 = jntvalues.at(1);
	double t3 = jntvalues.at(2);
	double t4 = jntvalues.at(3);
	double t5 = jntvalues.at(4);
	double t6 = jntvalues.at(5);

	double t8 = sin(t4);
	double t9 = sin(t6);
	double t10 = cos(t4);
	double t11 = cos(t5);
	double t12 = cos(t6);
	double t13 = cos(t3);
	double t14 = t8*t9;
	double t30 = t10*t11*t12;
	double t15 = t14-t30;
	double t16 = sin(t3);
	double t17 = t9*t10;
	double t18 = t8*t11*t12;
	double t19 = t17+t18;
	double t20 = cos(t1);
	double t21 = cos(t2);
	double t22 = sin(t5);
	double t23 = sin(t2);
	double t24 = t8*9.97E-2;
	double t39 = t10*t22*9.96E-2;
	double t25 = t24-t39;
	double t26 = t10*9.97E-2;
	double t27 = t8*t22*9.96E-2;
	double t28 = t26+t27;
	double t29 = sin(t1);
	double t31 = t15*t16;
	double t77 = t13*t19;
	double t32 = t31-t77;
	double t33 = t23*t32;
	double t34 = t13*t15;
	double t35 = t16*t19;
	double t36 = t34+t35;
	double t37 = t33-t21*t36;
	double t38 = t16*t28;
	double t40 = t13*t25;
	double t57 = t13*3.922E-1;
	double t41 = t38+t40-t57;
	double t42 = t21*t41;
	double t43 = t16*3.922E-1;
	double t44 = t13*t28;
	double t58 = t16*t25;
	double t45 = t43+t44-t58;
	double t46 = t23*t45;
	double t56 = t21*(1.7E1/4.0E1);
	double t47 = t42+t46-t56;
	double t48 = t11*9.96E-2;
	double t49 = t48+1.333E-1;
	double t50 = t11*t12*9.97E-2;
	double t51 = t8*t12;
	double t52 = t9*t10*t11;
	double t53 = t51+t52;
	double t54 = t10*t12;
	double t61 = t8*t9*t11;
	double t55 = t54-t61;
	double t59 = t20*t49;
	double t60 = t59-t29*t47;
	double t62 = t13*t55;
	double t79 = t16*t53;
	double t63 = t62-t79;
	double t64 = t23*t63;
	double t65 = t13*t53;
	double t66 = t16*t55;
	double t67 = t65+t66;
	double t68 = t21*t67;
	double t69 = t64+t68;
	double t70 = t20*t47;
	double t71 = t29*t49;
	double t72 = t70+t71;
	double t73 = t8*t9*t11*3.922E-1;
	double t74 = t5-t6;
	double t75 = t5+t6;
	double t76 = t8*3.922E3;
	double t78 = t12*t22;
	J.at(0,0) = t60*(t20*t37+t12*t22*t29)+t72*(t29*t37-t12*t20*t22);
	J.at(0,1) = t50-t9*t10*3.922E-1-t9*t22*9.96E-2-t8*t11*t12*3.922E-1-t9*t10*t13*(1.7E1/4.0E1)+t8*t9*t16*(1.7E1/4.0E1)-t8*t11*t12*t13*(1.7E1/4.0E1)-t10*t11*t12*t16*(1.7E1/4.0E1);
	J.at(0,2) = t50-t9*t10*3.922E-1-t9*t22*9.96E-2-t8*t11*t12*3.922E-1;
	J.at(0,3) = cos(t74)/2.0E4+cos(t75)*9.965E-2;
	J.at(0,4) = t12*(-9.96E-2);
	J.at(1,0) = -t60*(t20*t69+t9*t22*t29)-t72*(t29*t69-t9*t20*t22);
	J.at(1,1) = t73-t9*t11*9.97E-2-t10*t12*3.922E-1-t12*t22*9.96E-2-t10*t12*t13*(1.7E1/4.0E1)+t8*t12*t16*(1.7E1/4.0E1)+t8*t9*t11*t13*(1.7E1/4.0E1)+t9*t10*t11*t16*(1.7E1/4.0E1);
	J.at(1,2) = t73-t9*t11*9.97E-2-t10*t12*3.922E-1-t12*t22*9.96E-2;
	J.at(1,3) = sin(t74)/2.0E4-sin(t75)*9.965E-2;
	J.at(1,4) = t9*9.96E-2;
	J.at(2,0) = sin(t2+t3+t4-t5)*1.68E-2+cos(t2-t5)*(1.7E1/8.0E1)+cos(t2+t3+t5)*1.961E-1+cos(t2+t3-t5)*1.961E-1-sin(t2+t3+t4+t5)*1.165E-1+cos(t2+t5)*(1.7E1/8.0E1);
	J.at(2,1) = (t22*(t76+sin(t3+t4)*4.25E3-9.97E2))/1.0E4;
	J.at(2,2) = (t22*(t76-9.97E2))/1.0E4;
	J.at(2,3) = t22*(-9.97E-2);
	J.at(3,0) = -t21*t32-t23*t36;
	J.at(3,1) = t78;
	J.at(3,2) = t78;
	J.at(3,3) = t78;
	J.at(3,4) = -t9;
	J.at(4,0) = t21*t63-t23*t67;
	J.at(4,1) = -t9*t22;
	J.at(4,2) = -t9*t22;
	J.at(4,3) = -t9*t22;
	J.at(4,4) = -t12;
	J.at(5,0) = -t22*sin(t2+t3+t4);
	J.at(5,1) = t11;
	J.at(5,2) = t11;
	J.at(5,3) = t11;
	J.at(5,5) = 1.0;

	return J;

}

//Jacobian expressed on the base frame
//// Simplified from Matlab
arma::mat PC::get_Jsym_spatial(const arma::vec jntvalues) {

	arma::mat J(6,6);
	J.zeros(); //very important to initialize to zero!

	double t1 = jntvalues.at(0);
	double t2 = jntvalues.at(1);
	double t3 = jntvalues.at(2);
	double t4 = jntvalues.at(3);
	double t5 = jntvalues.at(4);
	// double t6 = jntvalues.at(5);	

	double t7 = cos(t1);
	double t8 = sin(t1);
	double t9 = cos(t2);
	double t10 = cos(t3);
	double t11 = sin(t3);
	double t12 = cos(t4);
	double t13 = sin(t2);
	double t14 = sin(t4);
	double t15 = sin(t5);
	double t16 = t2+t3+t4;
	double t17 = cos(t16);
	double t18 = t17*9.97E2;
	double t19 = t2+t3+t4+t5;
	double t20 = cos(t19);
	double t21 = t2+t3;
	double t22 = sin(t21);
	double t23 = t22*3.922E3;
	double t24 = t2+t3+t4-t5;
	double t25 = cos(t24);
	double t26 = t25*4.98E2;
	double t27 = cos(t5);
	double t28 = t13*4.25E3;
	double t30 = t20*4.98E2;
	double t29 = t18+t23+t26+t28-t30;
	double t31 = sin(t24);
	double t32 = t31*4.98E-2;
	double t33 = sin(t16);
	double t34 = t33*9.97E-2;
	double t35 = sin(t19);
	double t36 = cos(t21);
	double t37 = t1+t2+t3+t4;
	double t38 = -t1+t2+t3+t4;

	J.at(0,0) = t7*1.333E-1+t8*t9*(1.7E1/4.0E1)+t7*t27*9.96E-2+t8*t9*t10*3.922E-1-t8*t11*t13*3.922E-1-t8*t9*t11*t12*9.97E-2-t8*t9*t10*t14*9.97E-2-t8*t10*t12*t13*9.97E-2+t8*t11*t13*t14*9.97E-2+t8*t9*t10*t12*t15*9.96E-2-t8*t9*t11*t14*t15*9.96E-2-t8*t11*t12*t13*t15*9.96E-2-t8*t10*t13*t14*t15*9.96E-2;
	J.at(0,1) = (t7*t29)/1.0E4;
	J.at(0,2) = (t7*(t18-t20*4.98E2+t23+t26))/1.0E4;
	J.at(0,3) = (t7*(t18-t20*4.98E2+t26))/1.0E4;
	J.at(0,4) = t8*t15*(-9.96E-2)-t7*t9*t10*t12*t27*9.96E-2+t7*t9*t11*t14*t27*9.96E-2+t7*t11*t12*t13*t27*9.96E-2+t7*t10*t13*t14*t27*9.96E-2;
	J.at(1,0) = t8*1.333E-1-t7*t9*(1.7E1/4.0E1)+t8*t27*9.96E-2-t7*t9*t10*3.922E-1+t7*t11*t13*3.922E-1+t7*t9*t11*t12*9.97E-2+t7*t9*t10*t14*9.97E-2+t7*t10*t12*t13*9.97E-2-t7*t11*t13*t14*9.97E-2-t7*t9*t10*t12*t15*9.96E-2+t7*t9*t11*t14*t15*9.96E-2+t7*t11*t12*t13*t15*9.96E-2+t7*t10*t13*t14*t15*9.96E-2;
	J.at(1,1) = (t8*t29)/1.0E4;
	J.at(1,2) = (t8*(t18+t23+t26-t30))/1.0E4;
	J.at(1,3) = (t8*(t18+t26-t30))/1.0E4;
	J.at(1,4) = t7*t15*9.96E-2-t8*t9*t10*t12*t27*9.96E-2+t8*t9*t11*t14*t27*9.96E-2+t8*t11*t12*t13*t27*9.96E-2+t8*t10*t13*t14*t27*9.96E-2;
	J.at(2,1) = t9*(-1.7E1/4.0E1)+t32+t34-t35*4.98E-2-t36*3.922E-1;
	J.at(2,2) = t32+t34-t35*4.98E-2-t36*3.922E-1;
	J.at(2,3) = t32+t34-t35*4.98E-2;
	J.at(2,4) = -t32-t35*4.98E-2;
	J.at(3,1) = t8;
	J.at(3,2) = t8;
	J.at(3,3) = t8;
	J.at(3,4) = sin(t37)/2.0+sin(t38)/2.0;
	J.at(3,5) = t8*t27-t7*t15*t17;
	J.at(4,1) = -t7;
	J.at(4,2) = -t7;
	J.at(4,3) = -t7;
	J.at(4,4) = cos(t37)*(-1.0/2.0)+cos(t38)/2.0;
	J.at(4,5) = -t7*t27-t8*t15*t17;
	J.at(5,0) = 1.0;
	J.at(5,4) = -t17;
	J.at(5,5) = -t15*t33;

	return J;
}